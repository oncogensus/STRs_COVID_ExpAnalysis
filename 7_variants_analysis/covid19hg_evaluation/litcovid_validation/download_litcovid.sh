#!/usr/bin/env bash
# download_litcovid.sh
# Baixa o arquivo gene-anotado do LitCovid (formato PubTator) para data/.
# (~2.3 GB; eh baixado uma unica vez e ignorado pelo git).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
mkdir -p data
URL="https://ftp.ncbi.nlm.nih.gov/pub/lu/LitCovid/litcovid2pubtator.json.gz"
OUT="data/litcovid2pubtator.json.gz"

if [ -s "$OUT" ] && [ "$(stat -c%s "$OUT" 2>/dev/null || echo 0)" -gt 100000000 ]; then
  echo "Ja existe e parece valido: $OUT (pule)"
  exit 0
fi

echo "Baixando $URL -> $OUT"
wget -O "$OUT" "$URL"

# validacao pos-download
if [ ! -s "$OUT" ] || [ "$(stat -c%s "$OUT")" -lt 100000000 ]; then
  echo "ERRO: download falhou ou arquivo suspeito (<100MB). Removendo." >&2
  rm -f "$OUT"
  exit 1
fi
echo "Pronto ($OUT)."
