#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${1:-}"
MICRON_DIR="${2:-}"
SAMPLES_CSV="${3:-}"
N_PER_UNIT="${4:-4}"

if [[ -z "$WORK_DIR" || -z "$MICRON_DIR" || -z "$SAMPLES_CSV" ]]; then
  echo "Usage: $0 WORK_DIR MICRODB_DIR SAMPLES_CSV [N_PER_UNIT]" >&2
  exit 2
fi

if [[ ! -d "$WORK_DIR" ]]; then
  echo "Missing WORK_DIR: $WORK_DIR" >&2
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
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

matches_sample_path_or_command() {
  local sample="$1"
  local path="$2"
  local dir cmd_sh

  if [[ "$path" == *"/${sample}/"* || "$path" == *"/${sample}."* || "$path" == *"_${sample}_"* ]]; then
    return 0
  fi

  dir="$(dirname "$path")"
  while [[ "$dir" == "$WORK_DIR"* && "$dir" != "/" ]]; do
    cmd_sh="$dir/.command.sh"
    if [[ -f "$cmd_sh" ]] && grep -Fq -- "--fasta_paths=${sample}.fasta" "$cmd_sh"; then
      return 0
    fi
    if [[ "$dir" == "$WORK_DIR" ]]; then
      break
    fi
    dir="$(dirname "$dir")"
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

cat_stream() {
  local p="$1"
  if [[ "$p" == *.zst ]]; then
    zstd -dc "$p"
  else
    cat "$p"
  fi
}

present_in_fasta() {
  local id="$1"
  local fasta="$2"
  awk -v id="$id" '
    /^>/ {
      h=substr($0,2)
      sub(/[[:space:]].*$/, "", h)
      if (h == id) found=1
    }
    END { exit(found ? 0 : 1) }
  ' "$fasta"
}

present_in_pdb() {
  local pdbid="$1"
  grep -iq "^${pdbid}" "$MICRON_DIR/pdb70/pdb70_a3m.ffindex" 2>/dev/null || \
    grep -iq "^${pdbid}" "$MICRON_DIR/pdb100/pdb100_2021Mar03_a3m.ffindex" 2>/dev/null
}

mode_find_expr() {
  case "$1" in
    uniref90) printf "%s" "\\( -name 'uniref90_hits.sto.zst' -o -name 'uniref90_hits.sto' \\)" ;;
    mgnify) printf "%s" "\\( -name 'mgnify_hits.sto.zst' -o -name 'mgnify_hits.sto' \\)" ;;
    small_bfd) printf "%s" "\\( -name 'small_bfd_hits.sto.zst' -o -name 'small_bfd_hits.sto' \\)" ;;
    pdb) printf "%s" "\\( -name 'pdb_hits.hhr' -o -name 't000_.hhr' \\)" ;;
    *) return 1 ;;
  esac
}

status=0

for sample in "${SAMPLES[@]}"; do
  echo "[$sample]"
  for mode in uniref90 mgnify small_bfd pdb; do
    list="$tmpdir/${sample}.${mode}.files"
    expr="$(mode_find_expr "$mode")"
    : > "$list"
    while IFS=$'\t' read -r ts path; do
      [[ -n "$path" ]] || continue
      if matches_sample_path_or_command "$sample" "$path"; then
        printf '%s\t%s\n' "$ts" "$path"
      fi
    done < <(eval "find \"$WORK_DIR\" -type f $expr -printf '%T@\\t%p\\n'") \
      | sort -nr | cut -f2- > "$list"

    if [[ ! -s "$list" ]]; then
      echo "  [$mode] no matching files"
      continue
    fi

    : > "$tmpdir/${sample}.${mode}.units"
    while read -r path; do
      sample_unit_from_path "$sample" "$path"
    done < "$list" | awk '!seen[$0]++' > "$tmpdir/${sample}.${mode}.units"

    while read -r unit; do
      [[ -n "$unit" ]] || continue
      unit_files="$tmpdir/${sample}.${mode}.$(echo "$unit" | tr ':' '_').files"
      grep -F "$(printf '\t')" "$list" >/dev/null 2>&1 || true
      : > "$unit_files"
      while read -r path; do
        [[ "$(sample_unit_from_path "$sample" "$path")" == "$unit" ]] && printf '%s\n' "$path" >> "$unit_files"
      done < "$list"

      ids="$tmpdir/${sample}.${mode}.$(echo "$unit" | tr ':' '_').ids"
      case "$mode" in
        pdb)
          while read -r hhr; do
            [[ -n "$hhr" && -f "$hhr" ]] || continue
            awk '
              /^[[:space:]]*[0-9]+[[:space:]]+/ {
                hit=$2
                pdb=toupper(substr(hit,1,4))
                if (pdb ~ /^[0-9A-Z]{4}$/) print pdb
              }
            ' "$hhr"
          done < "$unit_files" | awk '!seen[$0]++' | awk -v n="$N_PER_UNIT" 'NR<=n' > "$ids"
          ;;
        *)
          while read -r sto; do
            [[ -n "$sto" && -f "$sto" ]] || continue
            cat_stream "$sto"
          done < "$unit_files" \
            | awk -v q="$sample" '/^[^#\/[:space:]][^[:space:]]*[[:space:]]/ {
                id=$1
                if (id ~ /^chain_[A-Za-z0-9]+$/) next
                if (id ~ ("^" q "(\\.|$)")) next
                print id
              }' \
            | awk '!seen[$0]++' \
            | awk -v n="$N_PER_UNIT" 'NR<=n' > "$ids"
          ;;
      esac

      missing=0
      while read -r id; do
        [[ -n "$id" ]] || continue
        case "$mode" in
          uniref90)
            present_in_fasta "$id" "$MICRON_DIR/uniref90/uniref90.fasta" || {
              echo "  [$mode][$unit] missing $id"
              missing=1
            }
            ;;
          mgnify)
            present_in_fasta "$id" "$MICRON_DIR/mgnify/mgy_clusters.fa" || {
              echo "  [$mode][$unit] missing $id"
              missing=1
            }
            ;;
          small_bfd)
            present_in_fasta "$id" "$MICRON_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta" || {
              echo "  [$mode][$unit] missing $id"
              missing=1
            }
            ;;
          pdb)
            present_in_pdb "$id" || {
              echo "  [$mode][$unit] missing $id"
              missing=1
            }
            ;;
        esac
      done < "$ids"

      if [[ "$missing" -eq 0 ]]; then
        echo "  [$mode][$unit] ok"
      else
        status=1
      fi
    done < "$tmpdir/${sample}.${mode}.units"
  done
done

exit "$status"
