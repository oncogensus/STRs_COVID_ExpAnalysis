#!/usr/bin/env bash
# 3_run_all.sh — Gera scripts bash para IGV.js de cada STR com outlier.
# Le STRs_ID unicos de str_samples_bams.tsv e cria um .sh por STR.
# Uso: bash 3_run_all.sh
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente."; echo "Gere os BEDs primeiro: qsub 1_generate_beds.pbs"; exit 1; }

strs=($(awk -F'\t' 'NR>1{print $2}' "$TSV" | sort -u))
[ ${#strs[@]} -eq 0 ] && { echo "Nenhum STR encontrado no TSV."; exit 1; }

mkdir -p scripts

PORT=8201
for s in "${strs[@]}"; do
  safe=$(echo "$s" | tr ':' '_')
  script="scripts/igv_${safe}.sh"
  cat > "$script" <<EOF
#!/usr/bin/env bash
# IGV.js para ${s}
# Uso: bash $script
cd "$BASE"
bash 2_igv_variant.sh "$s" $PORT
EOF
  chmod +x "$script"
  echo "Gerado: $script (porta $PORT) -> $s"
  PORT=$((PORT + 1))
done

echo
echo "============================================================"
echo "Scripts gerados em scripts/ para ${#strs[@]} STRs."
echo "Para rodar todos:"
echo "  for f in scripts/igv_*.sh; do bash \"\$f\" & done"
echo "============================================================"
