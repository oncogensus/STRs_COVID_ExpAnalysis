#!/usr/bin/env bash
# Sobe os 8 servidores IGV.js (um por variante) de uma vez.
# Cada variante fica em sua propria porta (8201-8208).
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$BASE"

for g in ROBO2 ANK3 CDH12 NKAIN2 SEMA6D KCNH1 KCNQ5 ST6GALNAC3; do
  bash "igv_per_variant/$g.sh" &
done

echo
echo "============================================================"
echo "Servidores IGV.js nas portas 8201-8208 (uma por variante)."
echo "No PC:  ssh -L 8201-8208:localhost:8201-8208 Carlos_Chagas"
echo "Abra no navegador: http://localhost:8201  (ROBO2), 8202 (ANK3), ..."
echo "Ctrl+C encerra todos."
echo "============================================================"

wait
