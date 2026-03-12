#!/usr/bin/env bash
set -euo pipefail

# Build truncated micron DBs from per-sample hit artifacts.
#
# For each sample ID, this script takes the most recent hit files under WORK_DIR,
# extracts top N hit IDs per DB, unions IDs across samples (capped to MAX_TOTAL),
# and truncates selected DB files in MICRON_DIR.
#
# Usage:
#   scripts/truncate_micron_from_sample_hits.sh [WORK_DIR] [SOURCE_DB_DIR] [MICRON_DIR] [SAMPLES_CSV] [N_PER_SAMPLE] [MAX_TOTAL] [LOG_FILE] [WINDOW_MIN] [COLABFOLD_N]

WORK_DIR="${1:-/srv/scratch/z5378336/inputs/test_inputs/work/Mar-11}"
SOURCE_DB_DIR="${2:-old}"
MICRON_DIR="${3:-micron}"
SAMPLES_CSV="${4:-/srv/scratch/z5378336/inputs/test_inputs/fulltest.csv}"
N_PER_SAMPLE="${5:-4}"

if [[ ! "$N_PER_SAMPLE" =~ ^[0-9]+$ || "$N_PER_SAMPLE" -lt 1 ]]; then
  echo "N_PER_SAMPLE must be a positive integer (got: $N_PER_SAMPLE)" >&2
  exit 1
fi

if [[ ! -d "$WORK_DIR" ]]; then
  echo "Missing WORK_DIR: $WORK_DIR" >&2
  exit 1
fi
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
WINDOW_MIN="${8:-70}"
COLABFOLD_N="${9:-$MAX_TOTAL}"
MMSEQS_NO_INDEX="${MMSEQS_NO_INDEX:-}"
# ColabFold local search with --db-load-mode 0 expects prefilter-compatible
# indices, i.e. createindex --index-subset 0.
MMSEQS_INDEX_PAR="${MMSEQS_INDEX_PAR:---index-subset 0}"

if [[ ! "$WINDOW_MIN" =~ ^[0-9]+$ || "$WINDOW_MIN" -lt 1 ]]; then
  echo "WINDOW_MIN must be a positive integer (got: $WINDOW_MIN)" >&2
  exit 1
fi
if [[ ! "$COLABFOLD_N" =~ ^[0-9]+$ || "$COLABFOLD_N" -lt 1 ]]; then
  echo "COLABFOLD_N must be a positive integer (got: $COLABFOLD_N)" >&2
  exit 1
fi

mkdir -p "$(dirname "$LOG_FILE")"
if [[ ! -f "$LOG_FILE" ]]; then
  : > "$LOG_FILE"
fi

# Mirror stdout/stderr to a persistent run log.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "==== truncate_micron_from_sample_hits $(date '+%F %T') ===="
echo "log_file=$LOG_FILE"
echo "window_min=$WINDOW_MIN"
echo "colabfold_n=$COLABFOLD_N"

# Required for ffindex_order, and for reading .sto.zst.
module load hhsuite/3.3.0 >/dev/null 2>&1 || true
module load zstd >/dev/null 2>&1 || true
module load mmseqs2 >/dev/null 2>&1 || true

# Match ColabFold setup_databases.sh behavior: merged index artifacts.
export MMSEQS_FORCE_MERGE=1

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

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

recent_sample_files() {
  local sample="$1"
  local mode="$2"
  local out_list="$3"
  local find_expr

  case "$mode" in
    uniref90)
      find_expr="\\( -name 'uniref90_hits.sto.zst' -o -name 'uniref90_hits.sto' \\)"
      ;;
    mgnify)
      find_expr="\\( -name 'mgnify_hits.sto.zst' -o -name 'mgnify_hits.sto' \\)"
      ;;
    small_bfd)
      find_expr="\\( -name 'small_bfd_hits.sto.zst' -o -name 'small_bfd_hits.sto' \\)"
      ;;
    pdb)
      find_expr="-name 'pdb_hits.hhr'"
      ;;
    *)
      echo "Unsupported mode: $mode" >&2
      return 1
      ;;
  esac

  eval "find \"$WORK_DIR\" -type f $find_expr -mmin -\"$WINDOW_MIN\" -printf '%T@\\t%p\\n'" \
    | awk -F'\t' -v s="$sample" 'index($2,"/" s "/") || index($2,"/" s ".") || index($2,"_" s "_")' \
    | sort -nr | cut -f2- > "$out_list"
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

  if [[ ! -s "$sto_list" ]]; then
    : > "$out_ids"
    return
  fi

  while read -r sto_path; do
    [[ -n "$sto_path" && -f "$sto_path" ]] || continue
    cat_stream "$sto_path"
  done < "$sto_list" \
    | awk -v q="$sample" '/^[^#\/[:space:]][^[:space:]]*[[:space:]]/ {
        id=$1
        # Drop obvious query/self rows (sample IDs and chain labels from multimer STOs).
        if (id ~ /^chain_[A-Za-z0-9]+$/) next
        if (id ~ ("^" q "(\\.|$)")) next
        print id
      }' \
    | awk '!seen[$0]++' \
    | awk -v n="$N_PER_SAMPLE" 'NR<=n' > "$out_ids"
}

extract_pdb_ids() {
  local hhr_list="$1"
  local out_ids="$2"
  if [[ ! -s "$hhr_list" ]]; then
    : > "$out_ids"
    return
  fi

  while read -r hhr; do
    [[ -n "$hhr" && -f "$hhr" ]] || continue
    awk '
      /^[[:space:]]*[0-9]+[[:space:]]+/ {
        hit=$2
        pdb=toupper(substr(hit,1,4))
        if (pdb ~ /^[0-9A-Z]{4}$/) print pdb
      }
    ' "$hhr"
  done < "$hhr_list" \
    | awk '!seen[$0]++' \
    | awk -v n="$N_PER_SAMPLE" 'NR<=n' > "$out_ids"
}

recent_sample_a3m_files() {
  local sample="$1"
  local out_list="$2"
  find "$WORK_DIR" -type f -name '*.a3m' -mmin -"$WINDOW_MIN" -printf '%T@\t%p\n' \
    | awk -F'\t' -v s="$sample" 'index($2,"/" s "/") || index($2,"/" s ".") || index($2,"_" s "_")' \
    | sort -nr | cut -f2- > "$out_list"
}

extract_a3m_ids() {
  local a3m_list="$1"
  local sample="$2"
  local out_ids="$3"

  if [[ ! -s "$a3m_list" ]]; then
    : > "$out_ids"
    return
  fi

  while read -r a3m; do
    [[ -n "$a3m" && -f "$a3m" ]] || continue
    awk -v q="$sample" '
      /^>/ {
        id=substr($0,2)
        sub(/[[:space:]].*$/, "", id)
        if (id ~ /^chain_[A-Za-z0-9]+$/) next
        if (id ~ ("^" q "(\\.|$)")) next
        print id
      }
    ' "$a3m"
  done < "$a3m_list" \
    | awk 'NF>0 && !seen[$0]++' \
    > "$out_ids"
}

echo "Samples (${#SAMPLES[@]}): ${SAMPLES[*]}"
echo "N_PER_SAMPLE=$N_PER_SAMPLE MAX_TOTAL=$MAX_TOTAL"
echo "WINDOW_MIN=$WINDOW_MIN"
echo "COLABFOLD_N=$COLABFOLD_N"
echo "MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-0}"
echo "MMSEQS_INDEX_PAR=${MMSEQS_INDEX_PAR:-<none>}"

: > "$tmpdir/msa_hit_tokens.raw"

for sample in "${SAMPLES[@]}"; do
  u_list="$tmpdir/uniref90.${sample}.files"
  m_list="$tmpdir/mgnify.${sample}.files"
  s_list="$tmpdir/small_bfd.${sample}.files"
  p_list="$tmpdir/pdb.${sample}.files"
  a3m_list="$tmpdir/a3m.${sample}.files"

  recent_sample_files "$sample" uniref90 "$u_list"
  recent_sample_files "$sample" mgnify "$m_list"
  recent_sample_files "$sample" small_bfd "$s_list"
  recent_sample_files "$sample" pdb "$p_list"
  recent_sample_a3m_files "$sample" "$a3m_list"

  echo "[$sample]"
  echo "  uniref90 files(last ${WINDOW_MIN}m): $(wc -l < "$u_list")"
  sed 's/^/    /' "$u_list" | head -n 5 || true
  echo "  mgnify files(last ${WINDOW_MIN}m): $(wc -l < "$m_list")"
  sed 's/^/    /' "$m_list" | head -n 5 || true
  echo "  small_bfd files(last ${WINDOW_MIN}m): $(wc -l < "$s_list")"
  sed 's/^/    /' "$s_list" | head -n 5 || true
  echo "  pdb files(last ${WINDOW_MIN}m): $(wc -l < "$p_list")"
  sed 's/^/    /' "$p_list" | head -n 5 || true
  echo "  a3m files(last ${WINDOW_MIN}m): $(wc -l < "$a3m_list")"
  sed 's/^/    /' "$a3m_list" | head -n 5 || true

  extract_sto_ids "$u_list" "$sample" "$tmpdir/uniref90.${sample}.ids"
  extract_sto_ids "$m_list" "$sample" "$tmpdir/mgnify.${sample}.ids"
  extract_sto_ids "$s_list" "$sample" "$tmpdir/small_bfd.${sample}.ids"
  extract_pdb_ids "$p_list" "$tmpdir/pdb.${sample}.ids"
  extract_a3m_ids "$a3m_list" "$sample" "$tmpdir/a3m.${sample}.ids"

  cat "$tmpdir/uniref90.${sample}.ids" "$tmpdir/mgnify.${sample}.ids" "$tmpdir/small_bfd.${sample}.ids" "$tmpdir/a3m.${sample}.ids" >> "$tmpdir/msa_hit_tokens.raw"
done

awk 'NF>0 && !seen[$0]++' "$tmpdir/msa_hit_tokens.raw" > "$tmpdir/msa_hit_tokens.uniq"
echo "MSA-derived token count: $(wc -l < "$tmpdir/msa_hit_tokens.uniq")"

aggregate_ids() {
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

  cat "${files[@]}" \
    | awk 'NF>0 && !seen[$0]++' \
    | awk -v n="$MAX_TOTAL" 'NR<=n' > "$out"
}

aggregate_ids "$tmpdir/uniref90.*.ids" "$tmpdir/uniref90.ids"
aggregate_ids "$tmpdir/mgnify.*.ids" "$tmpdir/mgnify.ids"
aggregate_ids "$tmpdir/small_bfd.*.ids" "$tmpdir/small_bfd.ids"
aggregate_ids "$tmpdir/pdb.*.ids" "$tmpdir/pdb.ids"

echo "Aggregated IDs:"
for f in "$tmpdir"/uniref90.ids "$tmpdir"/mgnify.ids "$tmpdir"/small_bfd.ids "$tmpdir"/pdb.ids; do
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
      /^>/ {k++; keep=(k<=n)}
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
      /^>/ {k++; keep=(k<=n)}
      {if(keep) print}
    ' "$src_fa" > "$out_fa"
  fi
}

# FASTA updates.
pre_uniref90_sig="$(file_sig "$MICRON_DIR/uniref90/uniref90.fasta")"
pre_mgnify_sig="$(file_sig "$MICRON_DIR/mgnify/mgy_clusters.fa")"
pre_smallbfd_sig="$(file_sig "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

trim_fasta_by_ids_or_topn "$SOURCE_DB_DIR/uniref90/uniref90.fasta" "$MICRON_DIR/uniref90/uniref90.fasta" "$tmpdir/uniref90.ids" "$MAX_TOTAL"
trim_fasta_by_ids_or_topn "$SOURCE_DB_DIR/mgnify/mgy_clusters.fa" "$MICRON_DIR/mgnify/mgy_clusters.fa" "$tmpdir/mgnify.ids" "$MAX_TOTAL"
trim_fasta_by_ids_or_topn "$SOURCE_DB_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$tmpdir/small_bfd.ids" "$MAX_TOTAL"

post_uniref90_sig="$(file_sig "$MICRON_DIR/uniref90/uniref90.fasta")"
post_mgnify_sig="$(file_sig "$MICRON_DIR/mgnify/mgy_clusters.fa")"
post_smallbfd_sig="$(file_sig "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

record_change "$MICRON_DIR/uniref90/uniref90.fasta" "$pre_uniref90_sig" "$post_uniref90_sig" "records=$(fasta_records "$MICRON_DIR/uniref90/uniref90.fasta")"
record_change "$MICRON_DIR/mgnify/mgy_clusters.fa" "$pre_mgnify_sig" "$post_mgnify_sig" "records=$(fasta_records "$MICRON_DIR/mgnify/mgy_clusters.fa")"
record_change "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" "$pre_smallbfd_sig" "$post_smallbfd_sig" "records=$(fasta_records "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta")"

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
done < "$tmpdir/pdb.ids"
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
    awk -F'\t' 'NR==FNR {keep[$1]=1; print $1; next} !($1 in keep) {print $1}' "$ids" "$src_index" \
      | awk 'NF>0 && !seen[$0]++' \
      | awk -v n="$keep_n" 'NR<=n' > "$ids.tmp"
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

if command -v mmseqs >/dev/null 2>&1; then
  subset_mmseqs_db \
    "$SOURCE_DB_DIR/colabfold_uniref30" \
    "$MICRON_DIR/colabfold_uniref30" \
    "uniref30_2302_db" \
    "$COLABFOLD_N" \
    "$tmpdir/msa_hit_tokens.uniq"

  subset_mmseqs_db \
    "$SOURCE_DB_DIR/colabfold_envdb" \
    "$MICRON_DIR/colabfold_envdb" \
    "colabfold_envdb_202108_db" \
    "$COLABFOLD_N" \
    "$tmpdir/msa_hit_tokens.uniq"
else
  echo "WARN: mmseqs not found; skipping colabfold DB trimming" >&2
fi

post_cf_u30_bytes="$(dir_bytes "$MICRON_DIR/colabfold_uniref30")"
post_cf_env_bytes="$(dir_bytes "$MICRON_DIR/colabfold_envdb")"
record_change "$MICRON_DIR/colabfold_uniref30" "$pre_cf_u30_bytes" "$post_cf_u30_bytes" "bytes=$post_cf_u30_bytes"
record_change "$MICRON_DIR/colabfold_envdb" "$pre_cf_env_bytes" "$post_cf_env_bytes" "bytes=$post_cf_env_bytes"

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

awk 'NF>0 && !seen[$0]++' "$tmpdir/required_pdb.ids" > "$tmpdir/required_pdb.uniq.ids"

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
echo "Target scale: top ${N_PER_SAMPLE}/sample, capped at MAX_TOTAL=${MAX_TOTAL}."
