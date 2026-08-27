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
if [ -s "$OUT" ]; then
  echo "Ja existe: $OUT (pule)"
else
  echo "Baixando $URL -> $OUT"
  wget -O "$OUT" "$URL"
fi
echo "Pronto."
