#!/usr/bin/env bash
# run_all.sh — Sobe IGV.js para todos os genes com outlier.
# Verifica se data/str_samples_bams.tsv existe; se nao, avisa.
# Le genes unicos do TSV e inicia igv_variant.sh para cada um.
# Portas: 8201, 8202, ... (ordem alfabeta dos genes).
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente."; echo "Gere os BEDs primeiro: qsub 0_generate_beds.pbs"; exit 1; }

genes=($(awk -F'\t' 'NR>1{print $1}' "$TSV" | sort -u))
[ ${#genes[@]} -eq 0 ] && { echo "Nenhum gene encontrado no TSV."; exit 1; }

PORT=8201
for g in "${genes[@]}"; do
  echo "Iniciando $g na porta $PORT ..."
  bash "$(dirname "${BASH_SOURCE[0]}")/igv_variant.sh" "$g" "$PORT" &
  PORT=$((PORT + 1))
done

echo
echo "============================================================"
echo "Servidores IGV.js: ${#genes[@]} genes (portas 8201-$((PORT-1)))."
echo "No PC:  ssh -L 8201-$((PORT-1)):localhost:8201-$((PORT-1)) Carlos_Chagas"
echo "Ctrl+C encerra todos."
echo "============================================================"

wait
