#!/usr/bin/env bash
# 3_run_all.sh — Generates bash scripts for IGV.js for each STR with an outlier.
# Reads unique STRs_ID from str_samples_bams.tsv and creates one .sh per STR.
# Usage: bash 3_run_all.sh
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
[ -f "$TSV" ] || { echo "ERROR: $TSV missing."; echo "Generate the BEDs first: qsub 1_generate_beds.pbs"; exit 1; }

strs=($(awk -F'\t' 'NR>1{print $2}' "$TSV" | sort -u))
[ ${#strs[@]} -eq 0 ] && { echo "No STR found in TSV."; exit 1; }

mkdir -p scripts

PORT=8201
for s in "${strs[@]}"; do
  safe=$(echo "$s" | tr ':' '_')
  script="scripts/igv_${safe}.sh"
  cat > "$script" <<EOF
#!/usr/bin/env bash
# IGV.js for ${s}
# Usage: bash $script
cd "$BASE"
bash 2_igv_variant.sh "$s" $PORT
EOF
  chmod +x "$script"
  echo "Generated: $script (port $PORT) -> $s"
  PORT=$((PORT + 1))
done

echo
echo "============================================================"
echo "Scripts generated in scripts/ for ${#strs[@]} STRs."
echo "To run all:"
echo "  for f in scripts/igv_*.sh; do bash \"\$f\" & done"
echo "============================================================"
