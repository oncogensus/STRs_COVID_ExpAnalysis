#!/bin/bash
# Quick check for FOXB1 sources
set -euo pipefail

echo "=== FOXB1 ==="

# ReMap storage casings
for variant in "FOXB1" "foxb1" "Foxb1"; do
  url="https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$variant/remap2022_${variant}_all_macs2_hg38_v1_0.bed.gz"
  code=$(curl -sS -o /dev/null -w "%{http_code}" -I "$url" --max-time 10 2>/dev/null || echo 000)
  echo "  $variant: $code"
done

# ReMap target page
code=$(curl -sS -o /dev/null -w "%{http_code}" -I "https://remap.univ-amu.fr/target_page/FOXB1:9606" --max-time 10 2>/dev/null || echo 000)
echo "  Target page: $code"

# ENCODE search
response=$(curl -sS "https://www.encodeproject.org/search/?type=Experiment&target.title=FOXB1&assembly=GRCh38&status=released&format=json" --max-time 15 2>/dev/null || echo '{"@graph":[]}')
count=$(echo "$response" | jq '.["@graph"] | length' 2>/dev/null || echo 0)
echo "  ENCODE experiments: $count"
if [ "$count" -gt 0 ]; then
  echo "$response" | jq -r '.["@graph"][] | "\(.accession) \(.target.title) \(.biosample_ontology.term_name // "N/A")"' 2>/dev/null | head -10
fi

# Cistrome
code=$(curl -sS -o /dev/null -w "%{http_code}" -I "http://cistrome.org/db/#FOXB1" --max-time 10 2>/dev/null || echo 000)
echo "  Cistrome: $code"