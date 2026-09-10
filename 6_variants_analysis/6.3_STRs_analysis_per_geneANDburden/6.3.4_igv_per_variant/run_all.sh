#!/usr/bin/env bash
# run_all.sh — Gera scripts bash para IGV.js de cada gene com outlier.
# Le genes unicos de str_samples_bams.tsv e cria um .sh por gene.
# Uso: bash run_all.sh
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente."; echo "Gere os BEDs primeiro: qsub 0_generate_beds.pbs"; exit 1; }

genes=($(awk -F'\t' 'NR>1{print $1}' "$TSV" | sort -u))
[ ${#genes[@]} -eq 0 ] && { echo "Nenhum gene encontrado no TSV."; exit 1; }

mkdir -p scripts

PORT=8201
for g in "${genes[@]}"; do
  script="scripts/igv_${g}.sh"
  cat > "$script" <<EOF
#!/usr/bin/env bash
# IGV.js para ${g}
# Uso: bash $script
cd "$BASE"
bash igv_variant.sh "$g" $PORT
EOF
  chmod +x "$script"
  echo "Gerado: $script (porta $PORT)"
  PORT=$((PORT + 1))
done

echo
echo "============================================================"
echo "Scripts gerados em scripts/ para ${#genes[@]} genes."
echo "Para rodar todos:"
echo "  for f in scripts/igv_*.sh; do bash \"\$f\" & done"
echo "============================================================"
