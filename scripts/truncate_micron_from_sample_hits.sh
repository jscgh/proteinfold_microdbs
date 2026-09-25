#!/usr/bin/env bash
set -euo pipefail

# Build truncated micron DBs from per-sample hit artifacts.
#
# For each sample ID, this script takes the most recent hit files under WORK_DIR,
# extracts top N hit IDs per DB, unions IDs across samples (capped to MAX_TOTAL),
# and truncates selected DB files in MICRON_DIR.
#
# Usage:
#   scripts/truncate_micron_from_sample_hits.sh [WORK_DIR] [SOURCE_DB_DIR] [MICRON_DIR] [SAMPLES_CSV] [N_PER_SAMPLE] [MAX_TOTAL] [LOG_FILE] [WINDOW_MIN] [COLABFOLD_N] [PARAMS_SOURCE_DIR] [BOLTZ_MOLS_KEEP]
#
# Set WINDOW_MIN=0 to search the full archived WORK_DIR instead of only recent files.

WORK_DIR="${1:-}"
SOURCE_DB_DIR="${2:-}"
MICRON_DIR="${3:-}"
SAMPLES_CSV="${4:-}"
N_PER_SAMPLE="${5:-4}"

if [[ ! "$N_PER_SAMPLE" =~ ^[0-9]+$ || "$N_PER_SAMPLE" -lt 1 ]]; then
  echo "N_PER_SAMPLE must be a positive integer (got: $N_PER_SAMPLE)" >&2
  exit 1
fi

if [[ -z "$WORK_DIR" || -z "$SOURCE_DB_DIR" || -z "$MICRON_DIR" || -z "$SAMPLES_CSV" ]]; then
  echo "Usage: $0 WORK_DIR SOURCE_DB_DIR MICRODB_DIR SAMPLES_CSV [N_PER_SAMPLE] [MAX_TOTAL] [LOG_FILE] [WINDOW_MIN] [COLABFOLD_N] [PARAMS_SOURCE_DIR] [BOLTZ_MOLS_KEEP]" >&2
  exit 2
fi

IFS=':' read -r -a WORK_DIRS <<< "$WORK_DIR"
if [[ ${#WORK_DIRS[@]} -eq 0 ]]; then
  echo "Missing WORK_DIR: $WORK_DIR" >&2
  exit 1
fi
for work_root in "${WORK_DIRS[@]}"; do
  if [[ ! -d "$work_root" ]]; then
    echo "Missing WORK_DIR root: $work_root" >&2
    exit 1
  fi
done
if [[ ! -d "$SOURCE_DB_DIR" ]]; then
  echo "Missing SOURCE_DB_DIR: $SOURCE_DB_DIR" >&2
  exit 1
fi
if [[ ! -d "$MICRON_DIR" ]]; then
  echo "Missing MICRON_DIR: $MICRON_DIR" >&2
  exit 1
fi
if [[ ! -f "$SAMPLES_CSV" ]]; then
  echo "Missing SAMPLES_CSV: $SAMPLES_CSV" >&2
  exit 1
fi

mapfile -t SAMPLES < <(awk -F, 'NR>1 && $1 != "" {print $1}' "$SAMPLES_CSV")
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "No sample IDs found in: $SAMPLES_CSV" >&2
  exit 1
fi

DEFAULT_MAX_TOTAL=$((N_PER_SAMPLE * ${#SAMPLES[@]}))
MAX_TOTAL="${6:-$DEFAULT_MAX_TOTAL}"
if [[ ! "$MAX_TOTAL" =~ ^[0-9]+$ || "$MAX_TOTAL" -lt 1 ]]; then
  echo "MAX_TOTAL must be a positive integer (got: $MAX_TOTAL)" >&2
  exit 1
fi

LOG_FILE="${7:-$MICRON_DIR/truncate_micron_from_sample_hits.log}"
WINDOW_MIN="${8:-0}"
COLABFOLD_N="${9:-$MAX_TOTAL}"
PARAMS_SOURCE_DIR="${10:-}"
BOLTZ_MOLS_KEEP="${11:-0}"
INCLUDE_PARAMS="${INCLUDE_PARAMS:-0}"
MMSEQS_NO_INDEX="${MMSEQS_NO_INDEX:-}"
NCBI_TAXDUMP_DIR="${NCBI_TAXDUMP_DIR:-${TAXDUMP_DIR:-}}"
# ColabFold local search with --db-load-mode 0 expects prefilter-compatible
# indices, i.e. createindex --index-subset 0.
MMSEQS_INDEX_PAR="${MMSEQS_INDEX_PAR:---index-subset 0}"

if [[ ! "$WINDOW_MIN" =~ ^[0-9]+$ || "$WINDOW_MIN" -lt 0 ]]; then
  echo "WINDOW_MIN must be a non-negative integer (got: $WINDOW_MIN)" >&2
  exit 1
fi
if [[ ! "$COLABFOLD_N" =~ ^[0-9]+$ || "$COLABFOLD_N" -lt 1 ]]; then
  echo "COLABFOLD_N must be a positive integer (got: $COLABFOLD_N)" >&2
  exit 1
fi
if [[ ! "$BOLTZ_MOLS_KEEP" =~ ^[0-9]+$ || "$BOLTZ_MOLS_KEEP" -lt 0 ]]; then
  echo "BOLTZ_MOLS_KEEP must be a non-negative integer (got: $BOLTZ_MOLS_KEEP)" >&2
  exit 1
fi
if [[ "$INCLUDE_PARAMS" != "0" && "$INCLUDE_PARAMS" != "1" ]]; then
  echo "INCLUDE_PARAMS must be 0 or 1 (got: $INCLUDE_PARAMS)" >&2
  exit 1
fi

if [[ -z "$PARAMS_SOURCE_DIR" ]]; then
  if [[ -d "$SOURCE_DB_DIR/params" ]]; then
    PARAMS_SOURCE_DIR="$SOURCE_DB_DIR/params"
  else
    PARAMS_SOURCE_DIR="$SOURCE_DB_DIR/params"
  fi
fi

mkdir -p "$(dirname "$LOG_FILE")"
if [[ ! -f "$LOG_FILE" ]]; then
  : > "$LOG_FILE"
fi

# Mirror stdout/stderr to a persistent run log.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==== truncate_micron_from_sample_hits $(date '+%F %T') ===="
echo "log_file=$LOG_FILE"
echo "work_dirs=${WORK_DIRS[*]}"
echo "window_min=$WINDOW_MIN"
echo "colabfold_n=$COLABFOLD_N"
echo "params_source_dir=$PARAMS_SOURCE_DIR"
echo "include_params=$INCLUDE_PARAMS"
echo "boltz_mols_keep=$BOLTZ_MOLS_KEEP"

# Required for ffindex_order, and for reading .sto.zst.
module load hhsuite/3.3.0 >/dev/null 2>&1 || true
module load zstd >/dev/null 2>&1 || true
module load mmseqs2 >/dev/null 2>&1 || true
module load apptainer >/dev/null 2>&1 || true

# Match ColabFold setup_databases.sh behavior: merged index artifacts.
export MMSEQS_FORCE_MERGE=1

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

declare -A PATH_MODE_CACHE

CHANGE_LINES=()

file_sig() {
  local p="$1"
  if [[ -f "$p" ]]; then
    sha256sum "$p" | awk '{print $1}'
  else
    echo "__MISSING__"
  fi
}

dir_sig_cif() {
  local d="$1"
  local files
  files="$(find "$d" -maxdepth 1 -type f -name '*.cif' | sort)"
  if [[ -z "$files" ]]; then
    echo "__EMPTY__"
    return
  fi
  # Hash paths+contents deterministically.
  echo "$files" | while read -r f; do
    [[ -n "$f" ]] || continue
    sha256sum "$f"
  done | sha256sum | awk '{print $1}'
}

fasta_records() {
  local p="$1"
  if [[ -f "$p" ]]; then
    awk '/^>/{n++} END{print n+0}' "$p"
  else
    echo 0
  fi
}

dir_bytes() {
  local p="$1"
  if [[ -e "$p" ]]; then
    du -s -B1 "$p" | awk '{print $1}'
  else
    echo 0
  fi
}

record_change() {
  local path="$1"
  local before_sig="$2"
  local after_sig="$3"
  local detail="$4"
  if [[ "$before_sig" != "$after_sig" ]]; then
    CHANGE_LINES+=("CHANGED: $path -> $detail")
  fi
}

matches_sample_path_or_command() {
  local sample="$1"
  local path="$2"
  local dir cmd_sh work_root

  if [[ "$path" == *"/${sample}/"* || "$path" == *"/${sample}."* || "$path" == *"_${sample}_"* ]]; then
    return 0
  fi

  for work_root in "${WORK_DIRS[@]}"; do
    if [[ "$path" != "$work_root"* ]]; then
      continue
    fi
    dir="$(dirname "$path")"
    while [[ "$dir" == "$work_root"* && "$dir" != "/" ]]; do
      cmd_sh="$dir/.command.sh"
      if [[ -f "$cmd_sh" ]] && grep -Fq -- "--fasta_paths=${sample}.fasta" "$cmd_sh"; then
        return 0
      fi
      if [[ "$dir" == "$work_root" ]]; then
        break
      fi
      dir="$(dirname "$dir")"
    done
  done

  return 1
}

sample_unit_from_path() {
  local sample="$1"
  local path="$2"
  local unit="$sample"

  if [[ "$path" =~ /${sample}/([A-Z0-9]+)/ ]]; then
    unit="${sample}:${BASH_REMATCH[1]}"
  elif [[ "$path" =~ /msas/([A-Z0-9])/ ]]; then
    unit="${sample}:${BASH_REMATCH[1]}"
  fi

  printf '%s\n' "$unit"
}

unit_file_tag() {
  local unit="$1"
  printf '%s\n' "${unit//:/__}"
}

nearest_command_sh() {
  local path="$1"
  local dir work_root cmd_sh

  for work_root in "${WORK_DIRS[@]}"; do
    if [[ "$path" != "$work_root"* ]]; then
      continue
    fi
    dir="$(dirname "$path")"
    while [[ "$dir" == "$work_root"* && "$dir" != "/" ]]; do
      cmd_sh="$dir/.command.sh"
      if [[ -f "$cmd_sh" ]]; then
        printf '%s\n' "$cmd_sh"
        return 0
      fi
      if [[ "$dir" == "$work_root" ]]; then
        break
      fi
      dir="$(dirname "$dir")"
    done
  done

  return 1
}

infer_mode_for_path() {
  local path="$1"
  local cached="${PATH_MODE_CACHE[$path]:-}"
  local cmd_sh=""

  if [[ -n "$cached" ]]; then
    printf '%s\n' "$cached"
    return 0
  fi

  if [[ "$path" == *"helixfold3"* ]]; then
    cached="helixfold3"
  elif [[ "$path" == *"rosettafold_all_atom"* || "$path" == *"rfaa"* ]]; then
    cached="rosettafold_all_atom"
  elif [[ "$path" == *"colabfold"* ]]; then
    cached="colabfold"
  elif [[ "$path" == *"boltz"* ]]; then
    cached="boltz"
  elif cmd_sh="$(nearest_command_sh "$path" 2>/dev/null)"; then
    if grep -Eqi 'NFCORE_PROTEINFOLD:HELIXFOLD3:|(^|[^A-Za-z])helixfold3([^A-Za-z]|$)' "$cmd_sh"; then
      cached="helixfold3"
    elif grep -Eqi 'NFCORE_PROTEINFOLD:ROSETTAFOLD_ALL_ATOM:|(^|[^A-Za-z])(rfaa|rosettafold_all_atom)([^A-Za-z]|$)' "$cmd_sh"; then
      cached="rosettafold_all_atom"
    elif grep -Eqi 'NFCORE_PROTEINFOLD:BOLTZ:|(^|[^A-Za-z])boltz([^A-Za-z]|$)|--af3-json' "$cmd_sh"; then
      cached="boltz"
    elif grep -Eqi 'NFCORE_PROTEINFOLD:COLABFOLD:|colabfold_batch|colabfold_search' "$cmd_sh"; then
      cached="colabfold"
    elif grep -Eqi 'NFCORE_PROTEINFOLD:ALPHAFOLD2:|(^|[^A-Za-z])alphafold([^A-Za-z]|$)' "$cmd_sh"; then
      cached="alphafold2"
    else
      cached="unknown"
    fi
  else
    cached="unknown"
  fi

  PATH_MODE_CACHE[$path]="$cached"
  printf '%s\n' "$cached"
}

first_line() {
  local p="$1"
  [[ -f "$p" ]] || return 1
  awk 'NF>0 {print; exit}' "$p"
}

build_unit_file_lists() {
  local sample="$1"
  local input_list="$2"
  local unit_dir="$3"

  rm -rf "$unit_dir"
  mkdir -p "$unit_dir"
  : > "$unit_dir/unit.order"

  while read -r path; do
    [[ -n "$path" && -f "$path" ]] || continue
    local unit tag
    unit="$(sample_unit_from_path "$sample" "$path")"
    tag="$(unit_file_tag "$unit")"
    printf '%s\n' "$path" >> "$unit_dir/${tag}.files"
    printf '%s\n' "$unit" >> "$unit_dir/unit.order"
  done < "$input_list"
}

copy_file_or_dir() {
  local src="$1"
  local dst="$2"
  rm -rf "$dst"
  mkdir -p "$(dirname "$dst")"
  cp -a "$src" "$dst"
}

copy_local_file_or_dir() {
  local src="$1"
  local dst="$2"
  rm -rf "$dst"
  mkdir -p "$(dirname "$dst")"
  cp -aL "$src" "$dst"
}

build_portable_alphafold_params() {
  local src_root="$1"
  local dst_root="$2"
  local canonical_dir="$dst_root/alphafold_params"

  local remove_dirs=(
    alphafold_params
    alphafold_params_2021-07-14
    alphafold_params_2022-01-19
    alphafold_params_2022-03-02
    alphafold_params_2022-12-06
    alphafold_params_colab_2021-10-27
    alphafold_params_colab_2022-03-02
    alphafold_params_colab_2022-12-06
  )

  local entry
  for entry in "${remove_dirs[@]}"; do
    rm -rf "$dst_root/$entry"
  done

  mkdir -p "$canonical_dir"

  local license_src=""
  for license_src in \
    "$src_root/alphafold_params_2022-12-06/LICENSE" \
    "$src_root/alphafold_params_2021-07-14/LICENSE" \
    "$src_root/alphafold_params_colab_2022-12-06/LICENSE"; do
    if [[ -f "$license_src" ]]; then
      copy_local_file_or_dir "$license_src" "$canonical_dir/LICENSE"
      break
    fi
  done

  local canonical_sources=(
    "params_model_1.npz:$src_root/alphafold_params_2022-12-06/params_model_1.npz"
    "params_model_1_ptm.npz:$src_root/alphafold_params_2022-12-06/params_model_1_ptm.npz"
    "params_model_1_multimer_v2.npz:$src_root/alphafold_params_colab_2022-12-06/params_model_1_multimer_v2.npz"
    "params_model_1_multimer_v3.npz:$src_root/alphafold_params_2022-12-06/params_model_1_multimer_v3.npz"
  )

  local pair dst_name src_path
  for pair in "${canonical_sources[@]}"; do
    dst_name="${pair%%:*}"
    src_path="${pair#*:}"
    if [[ -f "$src_path" ]]; then
      copy_local_file_or_dir "$src_path" "$canonical_dir/$dst_name"
    else
      echo "WARN: missing AlphaFold source file: $src_path" >&2
    fi
  done

  if [[ -f "$canonical_dir/params_model_1_multimer_v2.npz" ]]; then
    rm -f "$canonical_dir/params_model_1_multimer.npz"
    ln "$canonical_dir/params_model_1_multimer_v2.npz" \
      "$canonical_dir/params_model_1_multimer.npz"
  fi

  make_family_hardlinks() {
    local family_suffix="$1"
    local src_name="$canonical_dir/params_model_1${family_suffix}.npz"
    local model dst_name
    [[ -f "$src_name" ]] || return 0
    for model in 2 3 4 5; do
      dst_name="$canonical_dir/params_model_${model}${family_suffix}.npz"
      rm -f "$dst_name"
      ln "$src_name" "$dst_name"
    done
  }

  make_family_hardlinks ""
  make_family_hardlinks "_ptm"
  make_family_hardlinks "_multimer"
  make_family_hardlinks "_multimer_v2"
  make_family_hardlinks "_multimer_v3"

  local variant_link
  for variant_link in \
    alphafold_params_2021-07-14 \
    alphafold_params_2022-01-19 \
    alphafold_params_2022-03-02 \
    alphafold_params_2022-12-06 \
    alphafold_params_colab_2021-10-27 \
    alphafold_params_colab_2022-03-02 \
    alphafold_params_colab_2022-12-06; do
    rm -rf "$dst_root/$variant_link"
    ln -sfn "alphafold_params" "$dst_root/$variant_link"
  done
}

sync_params_for_boltz() {
  local src_root="$1"
  local dst_root="$2"

  if [[ ! -d "$src_root" ]]; then
    echo "WARN: missing PARAMS_SOURCE_DIR: $src_root (skipping params sync)" >&2
    return
  fi

  local pre_params_bytes post_params_bytes
  pre_params_bytes="$(dir_bytes "$dst_root")"

  if [[ -L "$dst_root" ]]; then
    echo "INFO: replacing params symlink with local directory at $dst_root" >&2
    rm -f "$dst_root"
  fi
  mkdir -p "$dst_root"

  local param_entries=(
    boltz1_conf.ckpt
    boltz2_aff.ckpt
    boltz2_conf.ckpt
    ccd.pkl
    ccd_preprocessed_etkdg.pkl.gz
    esm2_t36_3B_UR50D-contact-regression.pt
    esm2_t36_3B_UR50D.pt
    esmfold_3B_v1.pt
    HelixFold3-240814.pdparams
    RFAA_paper_weights.pt
  )
  # Boltz always loads its canonical CCD component set from params/mols.
  # Keep these even for "micro" databases so protein inputs do not fail with
  # errors like "CCD component ALA not found".
  local canonical_mol_entries=(
    ALA
    ARG
    ASN
    ASP
    CYS
    GLN
    GLU
    GLY
    HIS
    ILE
    LEU
    LYS
    MET
    PHE
    PRO
    SER
    THR
    TRP
    TYR
    VAL
    UNK
  )
  local remove_entries=(
    alphafold_params
    alphafold_params_2021-07-14
    alphafold_params_2022-01-19
    alphafold_params_2022-03-02
    alphafold_params_2022-12-06
    alphafold_params_colab_2021-10-27
    alphafold_params_colab_2022-03-02
    alphafold_params_colab_2022-12-06
  )

  local entry
  for entry in "${remove_entries[@]}"; do
    if [[ -e "$dst_root/$entry" ]]; then
      rm -rf "$dst_root/$entry" 2>/dev/null || \
        echo "WARN: could not remove stale params entry: $dst_root/$entry" >&2
    fi
  done
  for entry in "${param_entries[@]}"; do
    if [[ -e "$src_root/$entry" ]]; then
      copy_local_file_or_dir "$src_root/$entry" "$dst_root/$entry"
    fi
  done

  build_portable_alphafold_params "$src_root" "$dst_root"

  rm -rf "$dst_root/mols"
  mkdir -p "$dst_root/mols"
  if [[ -d "$src_root/mols" ]]; then
    local mol_name
    local extra_kept=0
    for mol_name in "${canonical_mol_entries[@]}"; do
      if [[ -f "$src_root/mols/${mol_name}.pkl" ]]; then
        cp -a "$src_root/mols/${mol_name}.pkl" "$dst_root/mols/${mol_name}.pkl"
      else
        echo "WARN: missing canonical Boltz mol file: $src_root/mols/${mol_name}.pkl" >&2
      fi
    done

    if [[ "$BOLTZ_MOLS_KEEP" -gt 0 ]]; then
      find "$src_root/mols" -maxdepth 1 -type f -name '*.pkl' | LC_ALL=C sort \
        | while read -r mol_file; do
            [[ -n "$mol_file" ]] || continue
            base_name="$(basename "$mol_file")"
            if [[ -e "$dst_root/mols/$base_name" ]]; then
              continue
            fi
            cp -a "$mol_file" "$dst_root/mols/$base_name"
            extra_kept=$((extra_kept + 1))
            if [[ "$extra_kept" -ge "$BOLTZ_MOLS_KEEP" ]]; then
              break
            fi
          done
    fi
  fi

  post_params_bytes="$(dir_bytes "$dst_root")"
  local mols_count
  mols_count="$(find "$dst_root/mols" -maxdepth 1 -type f -name '*.pkl' | wc -l)"
  record_change "$dst_root" "$pre_params_bytes" "$post_params_bytes" "bytes=$post_params_bytes mols=$mols_count"
}

recent_sample_files() {
  local sample="$1"
  local mode="$2"
  local out_list="$3"
  local find_expr
  local work_root

  case "$mode" in
    uniref90)
      find_expr="\\( -name 'uniref90_hits.sto.zst' -o -name 'uniref90_hits.sto' \\)"
      ;;
    mgnify)
      find_expr="\\( -name 'mgnify_hits.sto.zst' -o -name 'mgnify_hits.sto' \\)"
      ;;
    uniprot)
      find_expr="\\( -name 'uniprot_hits.sto.zst' -o -name 'uniprot_hits.sto' \\)"
      ;;
    small_bfd)
      find_expr="\\( -name 'small_bfd_hits.sto.zst' -o -name 'small_bfd_hits.sto' \\)"
      ;;
    pdb)
      find_expr="\\( -name 'pdb_hits.hhr' -o -name 't000_.hhr' \\)"
      ;;
    hf3_pdb_hits)
      find_expr="\\( -name 'pdb_hits.sto' -o -name 'pdb_hits.sto.zst' \\)"
      ;;
    hf3_features)
      find_expr="\\( -path '*/msas/*/features.pkl' \\)"
      ;;
    *)
      echo "Unsupported mode: $mode" >&2
      return 1
      ;;
  esac

  : > "$out_list"
  for work_root in "${WORK_DIRS[@]}"; do
    if [[ "$WINDOW_MIN" -eq 0 ]]; then
      while IFS=$'\t' read -r ts path; do
        [[ -n "$path" ]] || continue
        if matches_sample_path_or_command "$sample" "$path"; then
          printf '%s\t%s\n' "$ts" "$path"
        fi
      done < <(eval "find \"$work_root\" -type f $find_expr -printf '%T@\\t%p\\n'")
    else
      while IFS=$'\t' read -r ts path; do
        [[ -n "$path" ]] || continue
        if matches_sample_path_or_command "$sample" "$path"; then
          printf '%s\t%s\n' "$ts" "$path"
        fi
      done < <(eval "find \"$work_root\" -type f $find_expr -mmin -\"$WINDOW_MIN\" -printf '%T@\\t%p\\n'")
    fi
  done | sort -nr | cut -f2- > "$out_list"
}

cat_stream() {
  local p="$1"
  if [[ "$p" == *.zst ]]; then
    zstd -dc "$p"
  else
    cat "$p"
  fi
}

extract_sto_ids() {
  local sto_list="$1"
  local sample="$2"
  local out_ids="$3"
  local out_must="${out_ids}.must"

  : > "$out_must"
  if [[ ! -s "$sto_list" ]]; then
    : > "$out_ids"
    return
  fi

  local unit_dir="$tmpdir/units.sto.${sample}"
  local unit tag unit_files
  build_unit_file_lists "$sample" "$sto_list" "$unit_dir"
  : > "$out_ids.raw"
  declare -A must_mode_seen=()

  while read -r unit; do
    [[ -n "$unit" ]] || continue
    tag="$(unit_file_tag "$unit")"
    unit_files="$unit_dir/${tag}.files"
    while read -r sto_path; do
      [[ -n "$sto_path" && -f "$sto_path" ]] || continue
      local mode file_ids first_id
      mode="$(infer_mode_for_path "$sto_path")"
      file_ids="$tmpdir/sto.$(basename "$out_ids").$(basename "$sto_path").ids"
      cat_stream "$sto_path" \
      | awk -v q="$sample" '/^[^#\/[:space:]][^[:space:]]*[[:space:]]/ {
          id=$1
          # Drop obvious query/self rows (sample IDs and chain labels from multimer STOs).
          if (id ~ /^chain_[A-Za-z0-9]+$/) next
          if (id ~ ("^" q "(\\.|$)")) next
          print id
        }' \
      | awk 'NF>0 && !seen[$0]++' > "$file_ids"
      awk -v n="$N_PER_SAMPLE" 'NR<=n' "$file_ids" >> "$out_ids.raw"
      if [[ "$mode" != "unknown" && -z "${must_mode_seen[$mode]:-}" ]]; then
        first_id="$(first_line "$file_ids" || true)"
        if [[ -n "$first_id" ]]; then
          printf '%s\n' "$first_id" >> "$out_must"
          must_mode_seen[$mode]=1
        fi
      fi
      rm -f "$file_ids"
    done < "$unit_files"
  done < <(awk 'NF>0 && !seen[$0]++' "$unit_dir/unit.order")

  awk 'NF>0 && !seen[$0]++' "$out_ids.raw" > "$out_ids"
  awk 'NF>0 && !seen[$0]++' "$out_must" > "${out_must}.uniq"
  mv "${out_must}.uniq" "$out_must"
}

extract_pdb_ids() {
  local hhr_list="$1"
  local sample="$2"
  local out_ids="$3"
  local out_must="${out_ids}.must"

  : > "$out_must"
  if [[ ! -s "$hhr_list" ]]; then
    : > "$out_ids"
    return
  fi

  local unit_dir="$tmpdir/units.pdb.${sample}"
  local unit tag unit_files
  build_unit_file_lists "$sample" "$hhr_list" "$unit_dir"
  : > "$out_ids.raw"
  declare -A must_mode_seen=()

  while read -r unit; do
    [[ -n "$unit" ]] || continue
    tag="$(unit_file_tag "$unit")"
    unit_files="$unit_dir/${tag}.files"
    while read -r hhr; do
      [[ -n "$hhr" && -f "$hhr" ]] || continue
      local mode file_ids first_id
      mode="$(infer_mode_for_path "$hhr")"
      file_ids="$tmpdir/pdb.$(basename "$out_ids").$(basename "$hhr").ids"
      awk '
        /^[[:space:]]*[0-9]+[[:space:]]+/ {
          hit=$2
          pdb=toupper(substr(hit,1,4))
          if (pdb ~ /^[0-9A-Z]{4}$/) print pdb
        }
      ' "$hhr" | awk 'NF>0 && !seen[$0]++' > "$file_ids"
      awk -v n="$N_PER_SAMPLE" 'NR<=n' "$file_ids" >> "$out_ids.raw"
      if [[ "$mode" != "unknown" && -z "${must_mode_seen[$mode]:-}" ]]; then
        first_id="$(first_line "$file_ids" || true)"
        if [[ -n "$first_id" ]]; then
          printf '%s\n' "$first_id" >> "$out_must"
          must_mode_seen[$mode]=1
        fi
      fi
      rm -f "$file_ids"
    done < "$unit_files"
  done < <(awk 'NF>0 && !seen[$0]++' "$unit_dir/unit.order")

  awk 'NF>0 && !seen[$0]++' "$out_ids.raw" > "$out_ids"
  awk 'NF>0 && !seen[$0]++' "$out_must" > "${out_must}.uniq"
  mv "${out_must}.uniq" "$out_must"
}

extract_pdb_chain_ids() {
  local hhr_list="$1"
  local sample="$2"
  local out_ids="$3"
  local out_must="${out_ids}.must"

  : > "$out_must"
  if [[ ! -s "$hhr_list" ]]; then
    : > "$out_ids"
    return
  fi

  local unit_dir="$tmpdir/units.pdbchain.${sample}"
  local unit tag unit_files
  build_unit_file_lists "$sample" "$hhr_list" "$unit_dir"
  : > "$out_ids.raw"
  declare -A must_mode_seen=()

  while read -r unit; do
    [[ -n "$unit" ]] || continue
    tag="$(unit_file_tag "$unit")"
    unit_files="$unit_dir/${tag}.files"
    while read -r hhr; do
      [[ -n "$hhr" && -f "$hhr" ]] || continue
      local mode file_ids first_id
      mode="$(infer_mode_for_path "$hhr")"
      file_ids="$tmpdir/pdbchain.$(basename "$out_ids").$(basename "$hhr").ids"
      awk '
        /^[[:space:]]*[0-9]+[[:space:]]+/ {
          hit=$2
          if (hit ~ /^[0-9A-Za-z]{4}_[A-Za-z0-9]+$/) print hit
        }
      ' "$hhr" | awk 'NF>0 && !seen[$0]++' > "$file_ids"
      awk -v n="$N_PER_SAMPLE" 'NR<=n' "$file_ids" >> "$out_ids.raw"
      if [[ "$mode" != "unknown" && -z "${must_mode_seen[$mode]:-}" ]]; then
        first_id="$(first_line "$file_ids" || true)"
        if [[ -n "$first_id" ]]; then
          printf '%s\n' "$first_id" >> "$out_must"
          must_mode_seen[$mode]=1
        fi
      fi
      rm -f "$file_ids"
    done < "$unit_files"
  done < <(awk 'NF>0 && !seen[$0]++' "$unit_dir/unit.order")

  awk 'NF>0 && !seen[$0]++' "$out_ids.raw" > "$out_ids"
  awk 'NF>0 && !seen[$0]++' "$out_must" > "${out_must}.uniq"
  mv "${out_must}.uniq" "$out_must"
}

extract_hf3_pdb_hit_ids() {
  local sto_list="$1"
  local sample="$2"
  local out_full_ids="$3"
  local out_pdb_ids="$4"
  local out_full_must="${out_full_ids}.must"
  local out_pdb_must="${out_pdb_ids}.must"

  : > "$out_full_ids"
  : > "$out_pdb_ids"
  : > "$out_full_must"
  : > "$out_pdb_must"
  if [[ ! -s "$sto_list" ]]; then
    return
  fi

  local unit_dir="$tmpdir/units.hf3pdb.${sample}"
  local unit tag unit_files
  build_unit_file_lists "$sample" "$sto_list" "$unit_dir"
  : > "$out_full_ids.raw"
  : > "$out_pdb_ids.raw"

  while read -r unit; do
    [[ -n "$unit" ]] || continue
    tag="$(unit_file_tag "$unit")"
    unit_files="$unit_dir/${tag}.files"
    while read -r sto_path; do
      [[ -n "$sto_path" && -f "$sto_path" ]] || continue
      cat_stream "$sto_path"
    done < "$unit_files" \
      | awk '
          /^[#\/]/ { next }
          /^[[:space:]]*$/ { next }
          {
            id=$1
            if (id ~ /^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]+/) {
              print id > full_out
              pdb=toupper(substr(id,1,4))
              if (pdb ~ /^[0-9A-Z]{4}$/) print pdb > pdb_out
            }
          }
        ' full_out="$out_full_ids.raw" pdb_out="$out_pdb_ids.raw"
  done < <(awk 'NF>0 && !seen[$0]++' "$unit_dir/unit.order")

  awk 'NF>0 && !seen[$0]++' "$out_full_ids.raw" > "$out_full_ids"
  awk 'NF>0 && !seen[$0]++' "$out_pdb_ids.raw" > "$out_pdb_ids"
  first_line "$out_full_ids" > "$out_full_must" || true
  first_line "$out_pdb_ids" > "$out_pdb_must" || true
}

extract_hf3_selected_template_ids() {
  local features_list="$1"
  local sample="$2"
  local out_full_ids="$3"
  local out_pdb_ids="$4"
  local out_full_must="${out_full_ids}.must"
  local out_pdb_must="${out_pdb_ids}.must"

  : > "$out_full_ids"
  : > "$out_pdb_ids"
  : > "$out_full_must"
  : > "$out_pdb_must"
  if [[ ! -s "$features_list" ]]; then
    return
  fi

  local container_image="${HELIXFOLD3_CONTAINER_IMAGE:-}"
  local apptainer_bin=""
  if command -v apptainer >/dev/null 2>&1; then
    apptainer_bin="$(command -v apptainer)"
  elif command -v singularity >/dev/null 2>&1; then
    apptainer_bin="$(command -v singularity)"
  fi

  local unit_dir="$tmpdir/units.hf3feat.${sample}"
  local unit tag unit_files
  build_unit_file_lists "$sample" "$features_list" "$unit_dir"
  : > "$out_full_ids.raw"
  : > "$out_pdb_ids.raw"

  while read -r unit; do
    [[ -n "$unit" ]] || continue
    tag="$(unit_file_tag "$unit")"
    unit_files="$unit_dir/${tag}.files"
    if [[ -n "$apptainer_bin" && -f "$container_image" ]]; then
      while read -r feat_path; do
        [[ -n "$feat_path" && -f "$feat_path" ]] || continue
        "$apptainer_bin" exec "$container_image" python3.9 - <<'PY' "$feat_path"
import pickle, sys
path = sys.argv[1]
with open(path, 'rb') as handle:
    data = pickle.load(handle)
for raw_name in data.get('template_domain_names', []):
    if isinstance(raw_name, (bytes, bytearray)):
        name = raw_name.decode('utf-8', 'ignore')
    else:
        name = str(raw_name)
    if not name:
        continue
    print(name)
PY
      done < "$unit_files" \
        | awk 'NF>0 && !seen[$0]++' \
        | awk -v n="$N_PER_SAMPLE" 'NR<=n' >> "$out_full_ids.raw"
    fi
  done < <(awk 'NF>0 && !seen[$0]++' "$unit_dir/unit.order")

  if [[ ! -s "$out_full_ids.raw" ]]; then
    echo "WARN: failed to extract HF3 template_domain_names for $sample; falling back to pdb_hits.sto" >&2
    return 1
  fi

  awk 'NF>0 && !seen[$0]++' "$out_full_ids.raw" > "$out_full_ids"
  awk '
    NF>0 {
      pdb=toupper(substr($0,1,4))
      if (pdb ~ /^[0-9A-Z]{4}$/ && !seen[pdb]++) print pdb
    }
  ' "$out_full_ids" > "$out_pdb_ids"
  first_line "$out_full_ids" > "$out_full_must" || true
  first_line "$out_pdb_ids" > "$out_pdb_must" || true
}

recent_sample_a3m_files() {
  local sample="$1"
  local out_list="$2"
  local work_root
  : > "$out_list"
  for work_root in "${WORK_DIRS[@]}"; do
    if [[ "$WINDOW_MIN" -eq 0 ]]; then
      while IFS=$'\t' read -r ts path; do
        [[ -n "$path" ]] || continue
        if matches_sample_path_or_command "$sample" "$path"; then
          printf '%s\t%s\n' "$ts" "$path"
        fi
      done < <(find "$work_root" -type f -name '*.a3m' -printf '%T@\t%p\n')
    else
      while IFS=$'\t' read -r ts path; do
        [[ -n "$path" ]] || continue
        if matches_sample_path_or_command "$sample" "$path"; then
          printf '%s\t%s\n' "$ts" "$path"
        fi
      done < <(find "$work_root" -type f -name '*.a3m' -mmin "-$WINDOW_MIN" -printf '%T@\t%p\n')
    fi
  done | sort -nr | cut -f2- > "$out_list"
}

extract_a3m_ids() {
  local a3m_list="$1"
  local sample="$2"
  local out_ids="$3"
  local out_must="${out_ids}.must"

  : > "$out_must"
  if [[ ! -s "$a3m_list" ]]; then
    : > "$out_ids"
    return
  fi

  declare -A must_mode_seen=()
  while read -r a3m; do
    [[ -n "$a3m" && -f "$a3m" ]] || continue
    local mode file_ids first_id
    mode="$(infer_mode_for_path "$a3m")"
    file_ids="$tmpdir/a3m.$(basename "$out_ids").$(basename "$a3m").ids"
    awk -v q="$sample" '
      /^>/ {
        id=substr($0,2)
        sub(/[[:space:]].*$/, "", id)
        if (id ~ /^chain_[A-Za-z0-9]+$/) next
        if (id ~ ("^" q "(\\.|$)")) next
        print id
      }
    ' "$a3m" | awk 'NF>0 && !seen[$0]++' > "$file_ids"
    cat "$file_ids" >> "$out_ids.raw"
    if [[ "$mode" != "unknown" && -z "${must_mode_seen[$mode]:-}" ]]; then
      first_id="$(first_line "$file_ids" || true)"
      if [[ -n "$first_id" ]]; then
        printf '%s\n' "$first_id" >> "$out_must"
        must_mode_seen[$mode]=1
      fi
    fi
    rm -f "$file_ids"
  done < "$a3m_list"

  awk 'NF>0 && !seen[$0]++' "$out_ids.raw" > "$out_ids"
  awk 'NF>0 && !seen[$0]++' "$out_must" > "${out_must}.uniq"
  mv "${out_must}.uniq" "$out_must"
}

echo "Samples (${#SAMPLES[@]}): ${SAMPLES[*]}"
echo "N_PER_SAMPLE=$N_PER_SAMPLE MAX_TOTAL=$MAX_TOTAL"
echo "WINDOW_MIN=$WINDOW_MIN"
echo "COLABFOLD_N=$COLABFOLD_N"
echo "MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-0}"
echo "MMSEQS_INDEX_PAR=${MMSEQS_INDEX_PAR:-<none>}"
echo "NCBI_TAXDUMP_DIR=${NCBI_TAXDUMP_DIR:-<auto>}"
echo "PARAMS_SOURCE_DIR=$PARAMS_SOURCE_DIR"
echo "INCLUDE_PARAMS=$INCLUDE_PARAMS"
echo "BOLTZ_MOLS_KEEP=$BOLTZ_MOLS_KEEP"

: > "$tmpdir/msa_hit_tokens.raw"
: > "$tmpdir/sample_units.raw"

for sample in "${SAMPLES[@]}"; do
  u_list="$tmpdir/uniref90.${sample}.files"
  m_list="$tmpdir/mgnify.${sample}.files"
  up_list="$tmpdir/uniprot.${sample}.files"
  s_list="$tmpdir/small_bfd.${sample}.files"
  p_list="$tmpdir/pdb.${sample}.files"
  hf3_pdb_list="$tmpdir/hf3_pdb.${sample}.files"
  hf3_features_list="$tmpdir/hf3_features.${sample}.files"
  a3m_list="$tmpdir/a3m.${sample}.files"

  recent_sample_files "$sample" uniref90 "$u_list"
  recent_sample_files "$sample" mgnify "$m_list"
  recent_sample_files "$sample" uniprot "$up_list"
  recent_sample_files "$sample" small_bfd "$s_list"
  recent_sample_files "$sample" pdb "$p_list"
  recent_sample_files "$sample" hf3_pdb_hits "$hf3_pdb_list"
  recent_sample_files "$sample" hf3_features "$hf3_features_list"
  recent_sample_a3m_files "$sample" "$a3m_list"

  echo "[$sample]"
  echo "  uniref90 files(last ${WINDOW_MIN}m): $(wc -l < "$u_list")"
  sed 's/^/    /' "$u_list" | head -n 5 || true
  echo "  mgnify files(last ${WINDOW_MIN}m): $(wc -l < "$m_list")"
  sed 's/^/    /' "$m_list" | head -n 5 || true
  echo "  uniprot files(last ${WINDOW_MIN}m): $(wc -l < "$up_list")"
  sed 's/^/    /' "$up_list" | head -n 5 || true
  echo "  small_bfd files(last ${WINDOW_MIN}m): $(wc -l < "$s_list")"
  sed 's/^/    /' "$s_list" | head -n 5 || true
  echo "  pdb files(last ${WINDOW_MIN}m): $(wc -l < "$p_list")"
  sed 's/^/    /' "$p_list" | head -n 5 || true
  echo "  hf3 pdb_hits.sto files(last ${WINDOW_MIN}m): $(wc -l < "$hf3_pdb_list")"
  sed 's/^/    /' "$hf3_pdb_list" | head -n 5 || true
  echo "  hf3 features.pkl files(last ${WINDOW_MIN}m): $(wc -l < "$hf3_features_list")"
  sed 's/^/    /' "$hf3_features_list" | head -n 5 || true
  echo "  a3m files(last ${WINDOW_MIN}m): $(wc -l < "$a3m_list")"
  sed 's/^/    /' "$a3m_list" | head -n 5 || true

  extract_sto_ids "$u_list" "$sample" "$tmpdir/uniref90.${sample}.ids"
  extract_sto_ids "$m_list" "$sample" "$tmpdir/mgnify.${sample}.ids"
  extract_sto_ids "$up_list" "$sample" "$tmpdir/uniprot.${sample}.ids"
  extract_sto_ids "$s_list" "$sample" "$tmpdir/small_bfd.${sample}.ids"
  extract_pdb_ids "$p_list" "$sample" "$tmpdir/pdb.${sample}.ids"
  extract_pdb_chain_ids "$p_list" "$sample" "$tmpdir/pdb_chain.${sample}.ids"
  if ! extract_hf3_selected_template_ids "$hf3_features_list" "$sample" "$tmpdir/hf3_template_full.${sample}.ids" "$tmpdir/hf3_template_pdb.${sample}.ids"; then
    extract_hf3_pdb_hit_ids "$hf3_pdb_list" "$sample" "$tmpdir/hf3_template_full.${sample}.ids" "$tmpdir/hf3_template_pdb.${sample}.ids"
  fi
  extract_a3m_ids "$a3m_list" "$sample" "$tmpdir/a3m.${sample}.ids"

  cat "$u_list" "$m_list" "$up_list" "$s_list" "$p_list" "$hf3_pdb_list" "$hf3_features_list" "$a3m_list" 2>/dev/null \
    | while read -r path; do
        [[ -n "$path" && -f "$path" ]] || continue
        sample_unit_from_path "$sample" "$path"
      done >> "$tmpdir/sample_units.raw"

  cat "$tmpdir/uniref90.${sample}.ids" "$tmpdir/mgnify.${sample}.ids" "$tmpdir/small_bfd.${sample}.ids" "$tmpdir/a3m.${sample}.ids" >> "$tmpdir/msa_hit_tokens.raw"
done

awk 'NF>0 && !seen[$0]++' "$tmpdir/msa_hit_tokens.raw" > "$tmpdir/msa_hit_tokens.uniq"
echo "MSA-derived token count: $(wc -l < "$tmpdir/msa_hit_tokens.uniq")"
awk 'NF>0 && !seen[$0]++' "$tmpdir/sample_units.raw" > "$tmpdir/sample_units.uniq"
SAMPLE_UNIT_COUNT="$(wc -l < "$tmpdir/sample_units.uniq")"
REQUIRED_TOTAL=$((N_PER_SAMPLE * SAMPLE_UNIT_COUNT))
EFFECTIVE_MAX_TOTAL="$MAX_TOTAL"
if [[ "$REQUIRED_TOTAL" -gt "$EFFECTIVE_MAX_TOTAL" ]]; then
  EFFECTIVE_MAX_TOTAL="$REQUIRED_TOTAL"
fi
EFFECTIVE_COLABFOLD_N="$COLABFOLD_N"
if [[ "$REQUIRED_TOTAL" -gt "$EFFECTIVE_COLABFOLD_N" ]]; then
  EFFECTIVE_COLABFOLD_N="$REQUIRED_TOTAL"
fi
echo "Detected sample units ($SAMPLE_UNIT_COUNT): $(tr '\n' ' ' < "$tmpdir/sample_units.uniq")"
if [[ "$EFFECTIVE_MAX_TOTAL" != "$MAX_TOTAL" ]]; then
  echo "INFO: raising effective MAX_TOTAL from $MAX_TOTAL to $EFFECTIVE_MAX_TOTAL to preserve top ${N_PER_SAMPLE} hits across all discovered sample units"
fi
if [[ "$EFFECTIVE_COLABFOLD_N" != "$COLABFOLD_N" ]]; then
  echo "INFO: raising effective COLABFOLD_N from $COLABFOLD_N to $EFFECTIVE_COLABFOLD_N to preserve per-unit ColabFold coverage"
fi

aggregate_ids() {
  local glob="$1"
  local out="$2"

  local files=()
  local mandatory_files=()
  shopt -s nullglob
  files=( $glob )
  shopt -u nullglob

  if [[ ${#files[@]} -eq 0 ]]; then
    : > "$out"
    return
  fi

  local f
  for f in "${files[@]}"; do
    if [[ -s "${f}.must" ]]; then
      mandatory_files+=("${f}.must")
    fi
  done

  if [[ ${#mandatory_files[@]} -gt 0 ]]; then
    cat "${mandatory_files[@]}" | awk 'NF>0 && !seen[$0]++' > "$out.must"
  else
    : > "$out.must"
  fi

  cat "${files[@]}" | awk 'NF>0 && !seen[$0]++' > "$out.all"

  local must_count remaining
  must_count="$(wc -l < "$out.must" | tr -d ' ')"
  if [[ "$must_count" -ge "$EFFECTIVE_MAX_TOTAL" ]]; then
    cp -f "$out.must" "$out"
  else
    remaining=$((EFFECTIVE_MAX_TOTAL - must_count))
    cat "$out.must" <(awk 'NR==FNR {keep[$0]=1; next} !($0 in keep)' "$out.must" "$out.all" | awk -v n="$remaining" 'NR<=n') \
      | awk 'NF>0 && !seen[$0]++' > "$out"
  fi

  rm -f "$out.must" "$out.all"
}

aggregate_lists() {
  local glob="$1"
  local out="$2"

  local files=()
  shopt -s nullglob
  files=( $glob )
  shopt -u nullglob

  if [[ ${#files[@]} -eq 0 ]]; then
    : > "$out"
    return
  fi

  cat "${files[@]}" | awk 'NF>0 && !seen[$0]++' > "$out"
}

aggregate_ids "$tmpdir/uniref90.*.ids" "$tmpdir/uniref90.ids"
aggregate_ids "$tmpdir/mgnify.*.ids" "$tmpdir/mgnify.ids"
aggregate_ids "$tmpdir/uniprot.*.ids" "$tmpdir/uniprot.ids"
aggregate_ids "$tmpdir/small_bfd.*.ids" "$tmpdir/small_bfd.ids"
aggregate_ids "$tmpdir/pdb.*.ids" "$tmpdir/pdb.ids"
aggregate_ids "$tmpdir/pdb_chain.*.ids" "$tmpdir/pdb_chain.ids"
aggregate_ids "$tmpdir/hf3_template_pdb.*.ids" "$tmpdir/hf3_template_pdb.ids"
aggregate_ids "$tmpdir/hf3_template_full.*.ids" "$tmpdir/hf3_template_full.ids"
aggregate_lists "$tmpdir/uniref90.*.files" "$tmpdir/uniref90.files"
aggregate_lists "$tmpdir/mgnify.*.files" "$tmpdir/mgnify.files"
aggregate_lists "$tmpdir/uniprot.*.files" "$tmpdir/uniprot.files"
aggregate_lists "$tmpdir/small_bfd.*.files" "$tmpdir/small_bfd.files"

if [[ -s "$tmpdir/hf3_template_pdb.ids" ]]; then
  cat "$tmpdir/pdb.ids" "$tmpdir/hf3_template_pdb.ids" 2>/dev/null | awk 'NF>0 && !seen[$0]++' > "$tmpdir/pdb.all.ids"
else
  cp -f "$tmpdir/pdb.ids" "$tmpdir/pdb.all.ids"
fi

cat "$tmpdir/pdb_chain.ids" "$tmpdir/hf3_template_full.ids" 2>/dev/null | awk 'NF>0 && !seen[$0]++' > "$tmpdir/pdb_seqres.ids"

echo "Aggregated IDs:"
for f in "$tmpdir"/uniref90.ids "$tmpdir"/mgnify.ids "$tmpdir"/uniprot.ids "$tmpdir"/small_bfd.ids "$tmpdir"/pdb.ids "$tmpdir"/hf3_template_pdb.ids "$tmpdir"/hf3_template_full.ids "$tmpdir"/pdb_seqres.ids; do
  echo "  $(basename "$f"): $(wc -l < "$f")"
  sed 's/^/    /' "$f" || true
done

trim_fasta_by_ids_or_topn() {
  local src_fa="$1"
  local out_fa="$2"
  local ids_file="$3"
  local topn="$4"

  mkdir -p "$(dirname "$out_fa")"

  if [[ $(wc -l < "$ids_file") -eq 0 ]]; then
    awk -v n="$topn" '
      /^>/ {
        k++
        if (k > n) exit
        keep=1
      }
      {if(keep) print}
    ' "$src_fa" > "$out_fa"
    return
  fi

  # Read keep-IDs from ids_file first (NR==FNR) to avoid awk variable-file redirection issues.
  awk '
    NR==FNR {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      if ($0 != "") keepid[$0]=1
      next
    }
    /^>/ {
      hdr=substr($0,2)
      tok=hdr
      sub(/[[:space:]].*$/, "", tok)
      emit=(tok in keepid)
    }
    { if (emit) print }
  ' "$ids_file" "$src_fa" > "$out_fa"

  # Safety: if IDs existed but none matched headers, avoid empty output.
  if [[ $(awk '/^>/{n++} END{print n+0}' "$out_fa") -eq 0 ]]; then
    awk -v n="$topn" '
      /^>/ {
        k++
        if (k > n) exit
        keep=1
      }
      {if(keep) print}
    ' "$src_fa" > "$out_fa"
  fi
}

trim_fasta_from_sto_or_source() {
  local sto_list="$1"
  local src_fa="$2"
  local out_fa="$3"
  local ids_file="$4"
  local topn="$5"

  mkdir -p "$(dirname "$out_fa")"

  if [[ -s "$sto_list" && -s "$ids_file" ]]; then
    awk '
      NR==FNR {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
        if ($0 != "") {
          keep[$0]=1
          order[++n]=$0
        }
        next
      }
      /^[^#\/[:space:]][^[:space:]]*[[:space:]]/ {
        id=$1
        if (!(id in keep)) next
        seq=$2
        gsub(/[[:space:]]+/, "", seq)
        gsub(/[.-]/, "", seq)
        seqs[id]=seqs[id] seq
      }
      END {
        for (i=1; i<=n; i++) {
          id=order[i]
          if (id in seqs && seqs[id] != "") {
            print ">" id
            print seqs[id]
          }
        }
      }
    ' "$ids_file" <(while read -r sto_path; do
          [[ -n "$sto_path" && -f "$sto_path" ]] || continue
          cat_stream "$sto_path"
        done < "$sto_list") > "$out_fa"
  fi

  if [[ $(awk '/^>/{n++} END{print n+0}' "$out_fa" 2>/dev/null || echo 0) -eq 0 ]]; then
    trim_fasta_by_ids_or_topn "$src_fa" "$out_fa" "$ids_file" "$topn"
  fi
}

# FASTA updates.
pre_uniref90_sig="$(file_sig "$MICRON_DIR/uniref90/uniref90.fasta")"
pre_mgnify_sig="$(file_sig "$MICRON_DIR/mgnify/mgy_clusters.fa")"
pre_uniprot_sig="$(file_sig "$MICRON_DIR/uniprot/uniprot.fasta")"
pre_smallbfd_sig="$(file_sig "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

trim_fasta_from_sto_or_source "$tmpdir/uniref90.files" "$SOURCE_DB_DIR/uniref90/uniref90.fasta" "$MICRON_DIR/uniref90/uniref90.fasta" "$tmpdir/uniref90.ids" "$EFFECTIVE_MAX_TOTAL"
trim_fasta_from_sto_or_source "$tmpdir/mgnify.files" "$SOURCE_DB_DIR/mgnify/mgy_clusters.fa" "$MICRON_DIR/mgnify/mgy_clusters.fa" "$tmpdir/mgnify.ids" "$EFFECTIVE_MAX_TOTAL"
trim_fasta_from_sto_or_source "$tmpdir/uniprot.files" "$SOURCE_DB_DIR/uniprot/uniprot.fasta" "$MICRON_DIR/uniprot/uniprot.fasta" "$tmpdir/uniprot.ids" "$EFFECTIVE_MAX_TOTAL"
trim_fasta_from_sto_or_source "$tmpdir/small_bfd.files" "$SOURCE_DB_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$tmpdir/small_bfd.ids" "$EFFECTIVE_MAX_TOTAL"

post_uniref90_sig="$(file_sig "$MICRON_DIR/uniref90/uniref90.fasta")"
post_mgnify_sig="$(file_sig "$MICRON_DIR/mgnify/mgy_clusters.fa")"
post_uniprot_sig="$(file_sig "$MICRON_DIR/uniprot/uniprot.fasta")"
post_smallbfd_sig="$(file_sig "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

record_change "$MICRON_DIR/uniref90/uniref90.fasta" "$pre_uniref90_sig" "$post_uniref90_sig" "records=$(fasta_records "$MICRON_DIR/uniref90/uniref90.fasta")"
record_change "$MICRON_DIR/mgnify/mgy_clusters.fa" "$pre_mgnify_sig" "$post_mgnify_sig" "records=$(fasta_records "$MICRON_DIR/mgnify/mgy_clusters.fa")"
record_change "$MICRON_DIR/uniprot/uniprot.fasta" "$pre_uniprot_sig" "$post_uniprot_sig" "records=$(fasta_records "$MICRON_DIR/uniprot/uniprot.fasta")"
record_change "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$pre_smallbfd_sig" "$post_smallbfd_sig" "records=$(fasta_records "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

# HF3 template seqres FASTA.
mkdir -p "$MICRON_DIR/pdb_seqres"
pre_pdbseqres_sig="$(file_sig "$MICRON_DIR/pdb_seqres/pdb_seqres.txt")"
trim_fasta_by_ids_or_topn \
  "$SOURCE_DB_DIR/pdb_seqres/pdb_seqres.txt" \
  "$MICRON_DIR/pdb_seqres/pdb_seqres.txt" \
  "$tmpdir/pdb_seqres.ids" \
  "$EFFECTIVE_MAX_TOTAL"
post_pdbseqres_sig="$(file_sig "$MICRON_DIR/pdb_seqres/pdb_seqres.txt")"
record_change "$MICRON_DIR/pdb_seqres/pdb_seqres.txt" "$pre_pdbseqres_sig" "$post_pdbseqres_sig" "records=$(fasta_records "$MICRON_DIR/pdb_seqres/pdb_seqres.txt")"

# mmCIF files from aggregated pdb IDs.
mkdir -p "$MICRON_DIR/pdb_mmcif/mmcif_files"
pre_mmcif_sig="$(dir_sig_cif "$MICRON_DIR/pdb_mmcif/mmcif_files")"
while read -r pdbid; do
  [[ -n "$pdbid" ]] || continue
  lower_id="$(echo "$pdbid" | tr '[:upper:]' '[:lower:]')"
  src_cif="$SOURCE_DB_DIR/pdb_mmcif/mmcif_files/${lower_id}.cif"
  if [[ -f "$src_cif" ]]; then
    cp -a "$src_cif" "$MICRON_DIR/pdb_mmcif/mmcif_files/${lower_id}.cif"
  else
    echo "WARN: missing source CIF for $pdbid at $src_cif" >&2
  fi
done < "$tmpdir/pdb.all.ids"
mmcif_count="$(find "$MICRON_DIR/pdb_mmcif/mmcif_files" -maxdepth 1 -type f -name '*.cif' | wc -l)"

# pdb70 ffdb updates from aggregated pdb IDs.
if [[ -s "$tmpdir/pdb.ids" ]]; then
  keys="$tmpdir/pdb70.keys"
  : > "$keys"
  while read -r pdbid; do
    grep -i "^${pdbid}" "$SOURCE_DB_DIR/pdb70/pdb70_a3m.ffindex" | cut -f1 >> "$keys" || true
  done < "$tmpdir/pdb.ids"
  awk 'NF>0 && !seen[$0]++' "$keys" > "$tmpdir/pdb70.keys.uniq"

  if [[ -s "$tmpdir/pdb70.keys.uniq" ]]; then
    mkdir -p "$MICRON_DIR/pdb70"
    pre_pdb70_a3m_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_a3m.ffdata")"
    pre_pdb70_a3m_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_a3m.ffindex")"
    pre_pdb70_cs_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_cs219.ffdata")"
    pre_pdb70_cs_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_cs219.ffindex")"
    pre_pdb70_hhm_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_hhm.ffdata")"
    pre_pdb70_hhm_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_hhm.ffindex")"

    for t in a3m cs219 hhm; do
      ffindex_order \
        "$tmpdir/pdb70.keys.uniq" \
        "$SOURCE_DB_DIR/pdb70/pdb70_${t}.ffdata" \
        "$SOURCE_DB_DIR/pdb70/pdb70_${t}.ffindex" \
        "$MICRON_DIR/pdb70/pdb70_${t}.ffdata" \
        "$MICRON_DIR/pdb70/pdb70_${t}.ffindex"
    done

    post_pdb70_a3m_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_a3m.ffdata")"
    post_pdb70_a3m_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_a3m.ffindex")"
    post_pdb70_cs_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_cs219.ffdata")"
    post_pdb70_cs_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_cs219.ffindex")"
    post_pdb70_hhm_data_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_hhm.ffdata")"
    post_pdb70_hhm_idx_sig="$(file_sig "$MICRON_DIR/pdb70/pdb70_hhm.ffindex")"

    pdb70_entries="$(wc -l < "$MICRON_DIR/pdb70/pdb70_a3m.ffindex" 2>/dev/null || echo 0)"
    record_change "$MICRON_DIR/pdb70/pdb70_a3m.ffdata" "$pre_pdb70_a3m_data_sig" "$post_pdb70_a3m_data_sig" "ffindex_entries=$pdb70_entries"
    record_change "$MICRON_DIR/pdb70/pdb70_a3m.ffindex" "$pre_pdb70_a3m_idx_sig" "$post_pdb70_a3m_idx_sig" "ffindex_entries=$pdb70_entries"
    record_change "$MICRON_DIR/pdb70/pdb70_cs219.ffdata" "$pre_pdb70_cs_data_sig" "$post_pdb70_cs_data_sig" "ffindex_entries=$pdb70_entries"
    record_change "$MICRON_DIR/pdb70/pdb70_cs219.ffindex" "$pre_pdb70_cs_idx_sig" "$post_pdb70_cs_idx_sig" "ffindex_entries=$pdb70_entries"
    record_change "$MICRON_DIR/pdb70/pdb70_hhm.ffdata" "$pre_pdb70_hhm_data_sig" "$post_pdb70_hhm_data_sig" "ffindex_entries=$pdb70_entries"
    record_change "$MICRON_DIR/pdb70/pdb70_hhm.ffindex" "$pre_pdb70_hhm_idx_sig" "$post_pdb70_hhm_idx_sig" "ffindex_entries=$pdb70_entries"
  fi
fi

# pdb100 ffdb updates from aggregated pdb IDs.
if [[ -s "$tmpdir/pdb.ids" && -f "$SOURCE_DB_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex" ]]; then
  keys100="$tmpdir/pdb100.keys"
  : > "$keys100"
  while read -r pdbid; do
    grep -i "^${pdbid}" "$SOURCE_DB_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex" | cut -f1 >> "$keys100" || true
  done < "$tmpdir/pdb.ids"
  awk 'NF>0 && !seen[$0]++' "$keys100" > "$tmpdir/pdb100.keys.uniq"

  if [[ -s "$tmpdir/pdb100.keys.uniq" ]]; then
    mkdir -p "$MICRON_DIR/pdb100"
    pre_pdb100_a3m_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffdata")"
    pre_pdb100_a3m_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex")"
    pre_pdb100_cs_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffdata")"
    pre_pdb100_cs_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffindex")"
    pre_pdb100_hhm_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffdata")"
    pre_pdb100_hhm_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffindex")"
    pre_pdb100_pdb_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffdata")"
    pre_pdb100_pdb_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffindex")"

    for t in a3m cs219 hhm pdb; do
      ffindex_order \
        "$tmpdir/pdb100.keys.uniq" \
        "$SOURCE_DB_DIR/pdb100/pdb100_2021Mar03_${t}.ffdata" \
        "$SOURCE_DB_DIR/pdb100/pdb100_2021Mar03_${t}.ffindex" \
        "$MICRON_DIR/pdb100/pdb100_2021Mar03_${t}.ffdata" \
        "$MICRON_DIR/pdb100/pdb100_2021Mar03_${t}.ffindex"
    done

    post_pdb100_a3m_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffdata")"
    post_pdb100_a3m_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex")"
    post_pdb100_cs_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffdata")"
    post_pdb100_cs_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffindex")"
    post_pdb100_hhm_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffdata")"
    post_pdb100_hhm_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffindex")"
    post_pdb100_pdb_data_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffdata")"
    post_pdb100_pdb_idx_sig="$(file_sig "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffindex")"

    pdb100_entries="$(wc -l < "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex" 2>/dev/null || echo 0)"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffdata" "$pre_pdb100_a3m_data_sig" "$post_pdb100_a3m_data_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex" "$pre_pdb100_a3m_idx_sig" "$post_pdb100_a3m_idx_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffdata" "$pre_pdb100_cs_data_sig" "$post_pdb100_cs_data_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_cs219.ffindex" "$pre_pdb100_cs_idx_sig" "$post_pdb100_cs_idx_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffdata" "$pre_pdb100_hhm_data_sig" "$post_pdb100_hhm_data_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_hhm.ffindex" "$pre_pdb100_hhm_idx_sig" "$post_pdb100_hhm_idx_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffdata" "$pre_pdb100_pdb_data_sig" "$post_pdb100_pdb_data_sig" "ffindex_entries=$pdb100_entries"
    record_change "$MICRON_DIR/pdb100/pdb100_2021Mar03_pdb.ffindex" "$pre_pdb100_pdb_idx_sig" "$post_pdb100_pdb_idx_sig" "ffindex_entries=$pdb100_entries"
  fi
fi

# ColabFold DB trimming with mmseqs (optional but enabled by default).
subset_mmseqs_db() {
  local src_dir="$1"
  local dst_dir="$2"
  local db_prefix="$3"
  local keep_n="$4"
  local hits_file="$5"

  local src_index="$src_dir/${db_prefix}.index"
  if [[ ! -f "$src_index" ]]; then
    echo "WARN: missing mmseqs index: $src_index (skipping $db_prefix)" >&2
    return
  fi

  mkdir -p "$dst_dir"
  # Clear stale index artifacts so each run reflects current keep_n.
  rm -f \
    "$dst_dir/${db_prefix}.idx" \
    "$dst_dir/${db_prefix}.idx.dbtype" \
    "$dst_dir/${db_prefix}.idx.index" \
    "$dst_dir/${db_prefix}.idx.0" \
    "$dst_dir/${db_prefix}.idx.index.0"

  local ids="$tmpdir/${db_prefix}.ids"
  : > "$ids"

  # Prefer IDs mapped from previously identified MSA hits.
  local map_src=""
  if [[ -f "$src_dir/${db_prefix}.lookup" ]]; then
    map_src="$src_dir/${db_prefix}.lookup"
  elif [[ -f "$src_dir/${db_prefix%_db}_sample_seq.tsv" ]]; then
    map_src="$src_dir/${db_prefix%_db}_sample_seq.tsv"
  elif [[ -f "$src_dir/${db_prefix%_db}_sample.tsv" ]]; then
    map_src="$src_dir/${db_prefix%_db}_sample.tsv"
  fi

  if [[ -n "$map_src" && -s "$hits_file" ]]; then
    awk -F'\t' '
      NR==FNR {
        t=tolower($0)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
        if (t != "") hits[t]=1
        next
      }
      {
        id=$1
        hit_found=0
        for(i=2;i<=NF;i++) {
          f=tolower($i)
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", f)
          if (f in hits) {hit_found=1; break}
          split(f, a, /[|,; ]+/)
          for (k in a) {
            if (a[k] in hits) {hit_found=1; break}
          }
          if (hit_found) break
        }
        if (hit_found && id != "" && !seen[id]++) print id
      }
    ' "$hits_file" "$map_src" | awk -v n="$keep_n" 'NR<=n' > "$ids"
  fi

  # Fallback/supplement: top IDs from index.
  if [[ $(wc -l < "$ids") -lt "$keep_n" ]]; then
    cp "$ids" "$ids.tmp"
    awk -F'\t' -v n="$keep_n" 'NR<=n && $1 != "" {print $1}' "$src_index" >> "$ids.tmp"
    awk 'NF>0 && !seen[$0]++' "$ids.tmp" \
      | awk -v n="$keep_n" 'NR<=n' > "$ids.next"
    mv "$ids.next" "$ids.tmp"
    mv "$ids.tmp" "$ids"
  fi

  if [[ $(wc -l < "$ids") -eq 0 ]]; then
    echo "WARN: no ids selected from $src_index (skipping $db_prefix)" >&2
    return
  fi

  recreate_subdbs() {
    local suffix
    for suffix in "" "_h" "_seq" "_seq_h" "_aln"; do
      local sdb="$src_dir/${db_prefix}${suffix}"
      local ddb="$dst_dir/${db_prefix}${suffix}"
      # Remove prior output files so createsubdb cannot leave stale full DB pieces.
      rm -f "${ddb}" "${ddb}.index" "${ddb}.dbtype"
      rm -f "${ddb}."[0-9]*
      if [[ -f "${sdb}.index" && -f "${sdb}.dbtype" ]]; then
        mmseqs createsubdb "$ids" "$sdb" "$ddb" >/dev/null
      fi
    done
  }

  recreate_subdbs

  # Safety: ensure base DB cardinality is actually truncated to keep_n.
  local out_index="$dst_dir/${db_prefix}.index"
  if [[ -f "$out_index" ]]; then
    local out_n
    out_n="$(wc -l < "$out_index")"
    if [[ "$out_n" -gt "$keep_n" || "$out_n" -eq 0 ]]; then
      echo "WARN: $db_prefix subset cardinality unexpected (got $out_n, target <= $keep_n); forcing deterministic top-$keep_n rebuild" >&2
      awk -F'\t' -v n="$keep_n" 'NR<=n{print $1}' "$src_index" > "$ids"
      for suffix in "" "_h" "_seq" "_seq_h" "_aln"; do
        local sdb="$src_dir/${db_prefix}${suffix}"
        local ddb="$dst_dir/${db_prefix}${suffix}"
        rm -f "${ddb}" "${ddb}.index" "${ddb}.dbtype"
        rm -f "${ddb}."[0-9]*
        if [[ -f "${sdb}.index" && -f "${sdb}.dbtype" ]]; then
          mmseqs createsubdb "$ids" "$sdb" "$ddb" >/dev/null
        fi
      done
    fi
  fi

  trim_mapping_and_rebuild_taxonomy() {
    local mapping_src="$1"
    local lookup_src="$2"
    local mapping_dst="$dst_dir/${db_prefix}_mapping"
    local taxonomy_dst="$dst_dir/${db_prefix}_taxonomy"
    local trimmed_map="$tmpdir/${db_prefix}.mapping.trimmed"
    local taxdump_dir="${NCBI_TAXDUMP_DIR:-}"

    if [[ ! -f "$mapping_src" || ! -f "$lookup_src" ]]; then
      return
    fi

    awk '
      NR==FNR {
        keep[$1]=1
        next
      }
      ($1 in keep) { print }
    ' "$ids" "$mapping_src" > "$trimmed_map"

    if [[ -s "$trimmed_map" ]]; then
      cp -a "$trimmed_map" "$mapping_dst"
    fi

    if [[ -z "$taxdump_dir" ]]; then
      for candidate in \
        "$src_dir/taxonomy" \
        "$(dirname "$dst_dir")/taxonomy" \
        "$dst_dir/taxonomy" \
        "$SOURCE_DB_DIR/taxonomy" \
        "$SOURCE_DB_DIR/ncbi-taxdump"
      do
        if [[ -f "$candidate/nodes.dmp" && -f "$candidate/names.dmp" ]]; then
          taxdump_dir="$candidate"
          break
        fi
      done
    fi

    rm -f "$taxonomy_dst"
    if [[ -s "$mapping_dst" ]]; then
      if [[ -n "$taxdump_dir" && -f "$taxdump_dir/nodes.dmp" && -f "$taxdump_dir/names.dmp" ]]; then
        mmseqs createtaxdb \
          "$dst_dir/$db_prefix" \
          "$tmpdir/mmseqs_tax_${db_prefix}" \
          --ncbi-tax-dump "$taxdump_dir" \
          --tax-mapping-file "$mapping_dst" \
          --threads "${PBS_NCPUS:-8}" >/dev/null 2>&1 || \
          echo "WARN: createtaxdb failed for $db_prefix with taxdump at $taxdump_dir" >&2
      else
        mmseqs createtaxdb \
          "$dst_dir/$db_prefix" \
          "$tmpdir/mmseqs_tax_${db_prefix}" \
          --tax-mapping-file "$mapping_dst" \
          --threads "${PBS_NCPUS:-8}" >/dev/null 2>&1 || \
          echo "WARN: createtaxdb failed for $db_prefix without a local taxdump; kept trimmed mapping only" >&2
      fi
    fi
  }

  # Copy metadata/auxiliary files when present.
  for extra in \
    "${db_prefix}.lookup" \
    "${db_prefix}.GPU_READY" \
    "${db_prefix}_mapping" \
    "${db_prefix}_taxonomy" \
    "${db_prefix%_db}_sample.tsv" \
    "${db_prefix%_db}_sample_aln.tsv" \
    "${db_prefix%_db}_sample_h.tsv" \
    "${db_prefix%_db}_sample_seq.tsv"
  do
    if [[ -f "$src_dir/$extra" ]]; then
      # Skip self-copies when source/destination resolve to the same file.
      if [[ -e "$dst_dir/$extra" && "$src_dir/$extra" -ef "$dst_dir/$extra" ]]; then
        continue
      fi
      cp -a "$src_dir/$extra" "$dst_dir/$extra"
    fi
  done

  # If the source tree lacks these sidecars, preserve any pre-populated copies
  # in the destination tree and trim/rebuild them against the selected IDs.
  local mapping_candidate=""
  local lookup_candidate=""
  for candidate in \
    "$src_dir/${db_prefix}_mapping" \
    "$dst_dir/${db_prefix}_mapping"
  do
    if [[ -f "$candidate" ]]; then
      mapping_candidate="$candidate"
      break
    fi
  done
  for candidate in \
    "$src_dir/${db_prefix}.lookup" \
    "$dst_dir/${db_prefix}.lookup"
  do
    if [[ -f "$candidate" ]]; then
      lookup_candidate="$candidate"
      break
    fi
  done
  trim_mapping_and_rebuild_taxonomy "$mapping_candidate" "$lookup_candidate"

  if [[ -z "$MMSEQS_NO_INDEX" ]]; then
    # shellcheck disable=SC2086
    idx_log="$tmpdir/${db_prefix}.createindex.log"
    # shellcheck disable=SC2086
    mmseqs createindex "$dst_dir/$db_prefix" "$tmpdir/mmseqs_tmp_${db_prefix}" --remove-tmp-files 1 $MMSEQS_INDEX_PAR >"$idx_log" 2>&1 || true

    # For very tiny uniref30 subsets, prefilter kmers may be absent; grow until usable.
    if [[ "$db_prefix" == "uniref30_2302_db" ]] && grep -q "No k-mer could be extracted" "$idx_log"; then
      for grow_n in 64 128 256 512 1000; do
        if [[ "$grow_n" -le "$keep_n" ]]; then
          continue
        fi
        echo "WARN: $db_prefix lacks prefilter k-mers at keep_n=$keep_n; retrying with keep_n=$grow_n" >&2
        awk -F'\t' -v n="$grow_n" 'NR<=n{print $1}' "$src_index" > "$ids"
        recreate_subdbs
        # shellcheck disable=SC2086
        mmseqs createindex "$dst_dir/$db_prefix" "$tmpdir/mmseqs_tmp_${db_prefix}" --remove-tmp-files 1 $MMSEQS_INDEX_PAR >"$idx_log" 2>&1 || true
        if ! grep -q "No k-mer could be extracted" "$idx_log"; then
          keep_n="$grow_n"
          break
        fi
      done

      # Final safety fallback: keep ColabFold search functional by restoring the
      # full source uniref30 DB if subset indexing is still k-merless.
      if grep -q "No k-mer could be extracted" "$idx_log"; then
        echo "WARN: $db_prefix cannot build prefilter-capable index from subset; restoring full source DB for compatibility" >&2
        rm -f "$dst_dir/${db_prefix}"* "$dst_dir/${db_prefix}_"*
        cp -a "$src_dir/${db_prefix}"* "$dst_dir/"
      fi
    fi
  fi

  # ColabFold expects <db>.idx and <db>.idx.index. Some mmseqs builds emit
  # shard-suffixed files (<db>.idx.0, <db>.idx.index.0), so provide aliases.
  if [[ ! -e "$dst_dir/${db_prefix}.idx" && -e "$dst_dir/${db_prefix}.idx.0" ]]; then
    ln -s "${db_prefix}.idx.0" "$dst_dir/${db_prefix}.idx"
  fi
  if [[ ! -e "$dst_dir/${db_prefix}.idx.index" && -e "$dst_dir/${db_prefix}.idx.index.0" ]]; then
    ln -s "${db_prefix}.idx.index.0" "$dst_dir/${db_prefix}.idx.index"
  fi
  if [[ ! -e "$dst_dir/${db_prefix}.idx.dbtype" ]]; then
    if [[ -e "$src_dir/${db_prefix}.idx.dbtype" ]]; then
      cp -a "$src_dir/${db_prefix}.idx.dbtype" "$dst_dir/${db_prefix}.idx.dbtype"
    elif [[ -e "$dst_dir/${db_prefix}.dbtype" ]]; then
      cp -a "$dst_dir/${db_prefix}.dbtype" "$dst_dir/${db_prefix}.idx.dbtype"
    fi
  fi

  # Keep mapping/taxonomy aliases consistent with ColabFold setup output.
  if [[ -e "$dst_dir/${db_prefix}_mapping" ]]; then
    ln -sf "${db_prefix}_mapping" "$dst_dir/${db_prefix}.idx_mapping"
  fi
  if [[ -e "$dst_dir/${db_prefix}_taxonomy" ]]; then
    ln -sf "${db_prefix}_taxonomy" "$dst_dir/${db_prefix}.idx_taxonomy"
  fi
}

pre_cf_u30_bytes="$(dir_bytes "$MICRON_DIR/colabfold_uniref30")"
pre_cf_env_bytes="$(dir_bytes "$MICRON_DIR/colabfold_envdb")"

resolve_colabfold_source_dir() {
  local primary="$1"
  local fallback="$2"
  local db_prefix="$3"

  if [[ -f "$primary/${db_prefix}.index" ]]; then
    printf '%s\n' "$primary"
    return 0
  fi
  if [[ -f "$fallback/${db_prefix}.index" ]]; then
    printf '%s\n' "$fallback"
    return 0
  fi
  printf '%s\n' "$primary"
}

if command -v mmseqs >/dev/null 2>&1; then
  cf_u30_src="$(resolve_colabfold_source_dir "$SOURCE_DB_DIR/colabfold_uniref30" "$SOURCE_DB_DIR/colabfold_uniref30" "uniref30_2302_db")"
  cf_env_src="$(resolve_colabfold_source_dir "$SOURCE_DB_DIR/colabfold_envdb" "$SOURCE_DB_DIR/colabfold_envdb" "colabfold_envdb_202108_db")"

  subset_mmseqs_db \
    "$cf_u30_src" \
    "$MICRON_DIR/colabfold_uniref30" \
    "uniref30_2302_db" \
    "$EFFECTIVE_COLABFOLD_N" \
    "$tmpdir/msa_hit_tokens.uniq"

  subset_mmseqs_db \
    "$cf_env_src" \
    "$MICRON_DIR/colabfold_envdb" \
    "colabfold_envdb_202108_db" \
    "$EFFECTIVE_COLABFOLD_N" \
    "$tmpdir/msa_hit_tokens.uniq"
else
  echo "WARN: mmseqs not found; skipping colabfold DB trimming" >&2
fi

post_cf_u30_bytes="$(dir_bytes "$MICRON_DIR/colabfold_uniref30")"
post_cf_env_bytes="$(dir_bytes "$MICRON_DIR/colabfold_envdb")"
record_change "$MICRON_DIR/colabfold_uniref30" "$pre_cf_u30_bytes" "$post_cf_u30_bytes" "bytes=$post_cf_u30_bytes"
record_change "$MICRON_DIR/colabfold_envdb" "$pre_cf_env_bytes" "$post_cf_env_bytes" "bytes=$post_cf_env_bytes"

# Parameters are excluded by default because model weights are too large for
# this microdatabase repository and may have separate licence terms.
if [[ "$INCLUDE_PARAMS" == "1" ]]; then
  sync_params_for_boltz "$PARAMS_SOURCE_DIR" "$MICRON_DIR/params"
else
  echo "Skipping model parameters (set INCLUDE_PARAMS=1 to create a local, Git-ignored params tree)."
fi

# Keep mmCIF files aligned with what truncated pdb70/pdb100 indices actually reference.
: > "$tmpdir/required_pdb.ids"
for idx in \
  "$MICRON_DIR/pdb70/pdb70_a3m.ffindex" \
  "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex"
do
  if [[ -s "$idx" ]]; then
    awk -F'\t' '
      NF>0 {
        key=$1
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
        pdb=toupper(substr(key,1,4))
        if (pdb ~ /^[0-9A-Z]{4}$/) print pdb
      }
    ' "$idx" >> "$tmpdir/required_pdb.ids"
  fi
done

cat "$tmpdir/required_pdb.ids" "$tmpdir/hf3_template_pdb.ids" 2>/dev/null | awk 'NF>0 && !seen[$0]++' > "$tmpdir/required_pdb.merged.ids"
cp -f "$tmpdir/required_pdb.merged.ids" "$tmpdir/required_pdb.uniq.ids"

if [[ -s "$tmpdir/required_pdb.uniq.ids" ]]; then
  while read -r pdbid; do
    [[ -n "$pdbid" ]] || continue
    lower_id="$(echo "$pdbid" | tr '[:upper:]' '[:lower:]')"
    dst_cif="$MICRON_DIR/pdb_mmcif/mmcif_files/${lower_id}.cif"
    src_cif="$SOURCE_DB_DIR/pdb_mmcif/mmcif_files/${lower_id}.cif"
    if [[ ! -f "$dst_cif" ]]; then
      if [[ -f "$src_cif" ]]; then
        cp -a "$src_cif" "$dst_cif"
      else
        echo "WARN: template index references $pdbid but source CIF is missing: $src_cif" >&2
      fi
    fi
  done < "$tmpdir/required_pdb.uniq.ids"

  # Hard safety check: fail if any required CIF is still missing.
  missing_required=0
  while read -r pdbid; do
    [[ -n "$pdbid" ]] || continue
    lower_id="$(echo "$pdbid" | tr '[:upper:]' '[:lower:]')"
    dst_cif="$MICRON_DIR/pdb_mmcif/mmcif_files/${lower_id}.cif"
    if [[ ! -f "$dst_cif" ]]; then
      echo "ERROR: missing required CIF referenced by template index: $dst_cif" >&2
      missing_required=1
    fi
  done < "$tmpdir/required_pdb.uniq.ids"
  if [[ "$missing_required" -ne 0 ]]; then
    exit 1
  fi

  # Prune unreferenced mmCIF files (patch_mmcif.py-style behavior):
  # move extras to pdb_mmcif/orphaned for review.
  awk '{print tolower($0)}' "$tmpdir/required_pdb.uniq.ids" > "$tmpdir/required_pdb.lower.ids"
  orphan_dir="$MICRON_DIR/pdb_mmcif/orphaned"
  mkdir -p "$orphan_dir"
  moved_orphans=0
  while read -r cif_path; do
    [[ -n "$cif_path" ]] || continue
    cif_base="$(basename "$cif_path" .cif)"
    if ! grep -Fxq "$cif_base" "$tmpdir/required_pdb.lower.ids"; then
      mv "$cif_path" "$orphan_dir/$(basename "$cif_path")"
      moved_orphans=$((moved_orphans + 1))
    fi
  done < <(find "$MICRON_DIR/pdb_mmcif/mmcif_files" -maxdepth 1 -type f -name '*.cif' | sort)
  echo "Pruned mmCIF orphans moved: $moved_orphans -> $orphan_dir"
fi

post_mmcif_sig="$(dir_sig_cif "$MICRON_DIR/pdb_mmcif/mmcif_files")"
mmcif_count="$(find "$MICRON_DIR/pdb_mmcif/mmcif_files" -maxdepth 1 -type f -name '*.cif' | wc -l)"
record_change "$MICRON_DIR/pdb_mmcif/mmcif_files" "$pre_mmcif_sig" "$post_mmcif_sig" "cif_files=$mmcif_count"

echo "Change summary:"
if [[ ${#CHANGE_LINES[@]} -eq 0 ]]; then
  echo "  No output files changed."
else
  for line in "${CHANGE_LINES[@]}"; do
    echo "  $line"
  done
fi

echo "Done. Truncated DBs written to: $MICRON_DIR"
echo "Target scale: top ${N_PER_SAMPLE}/sample-unit, capped at MAX_TOTAL=${MAX_TOTAL} (effective ${EFFECTIVE_MAX_TOTAL})."
