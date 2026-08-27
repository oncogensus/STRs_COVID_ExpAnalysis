#!/usr/bin/env bash
# run_litcovid.sh  (Parte 1)
# Genes outlier do DBSCAN global -> literatura LitCovid.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
ROOT="$(cd "$HERE/../../.." && pwd)"
mkdir -p data results
CATALOG="${CATALOG:-$ROOT/samples/STRs_analysis_dataset.tsv}"
export CATALOG

echo "=== 0/2 baixa arquivo gene-anotado do LitCovid (se ausente) ==="
bash download_litcovid.sh

echo "=== 1/2 genes outlier do DBSCAN global (n_outliers>0) ==="
python3 get_outlier_genes.py --catalog "$CATALOG" --out data/outlier_genes.txt

echo "=== 2/2 literatura LitCovid (PubTator gene-anotado) ==="
python3 gene_literature.py --genes data/outlier_genes.txt \
  --litcovid data/litcovid2pubtator.json.gz \
  --out results/gene_literature.tsv \
  --summary-out results/gene_literature_summary.tsv

echo "Pronto. Veja results/gene_literature.tsv e results/gene_literature_summary.tsv"
