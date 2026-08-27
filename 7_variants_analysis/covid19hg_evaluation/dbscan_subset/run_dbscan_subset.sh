#!/usr/bin/env bash
# run_dbscan_subset.sh  (Parte 2)
# Genes sugestivos COVID-19 HG (p<1e-5) -> STRs na coorte -> re-roda DBSCAN no subset.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
ROOT="$(cd "$HERE/../../.." && pwd)"
mkdir -p data results
CATALOG="${CATALOG:-$ROOT/samples/STRs_analysis_dataset.tsv}"
RESIDUALS="${RESIDUALS:-$ROOT/5_dbscan/norm_test/STRs_normalized_residuals.tsv}"
WINDOW=50000
P_THRESH=1e-5

echo "=== 1/5 extrai genes sugestivos COVID-19 HG (p<$P_THRESH) ==="
python3 ../extract_covid_genes.py \
  --sumstats \
    ../data/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz \
  --gene-bed ../genes.hg38.bed \
  --out results/covid_genes_suggestive.tsv \
  --window "$WINDOW" --p-thresh "$P_THRESH"

echo "=== 2/5 overlap genes sugestivos x STRs da coorte ==="
python3 overlap_suggestive_strs.py \
  --catalog "$CATALOG" \
  --covid-genes results/covid_genes_suggestive.tsv \
  --out results/suggestive_gene_strs.tsv \
  --ids-out data/suggestive_strs_ids.txt \
  --window "$WINDOW"

echo "=== 3/5 subset de residuos + DBSCAN ==="
awk -F'\t' '
NR==FNR { ids[$1]; next }
FNR==1  { for(i=1;i<=NF;i++) if($i=="STRs_ID") c=i; print; next }
$(c) in ids
' data/suggestive_strs_ids.txt "$RESIDUALS" > data/suggestive_strs_residuals.tsv
Rscript run_dbscan_subset.R data/suggestive_strs_residuals.tsv results/suggestive_strs_outliers.tsv

echo "=== 4/5 junta outliers ==="
python3 summarize_subset.py \
  --overlap results/suggestive_gene_strs.tsv \
  --dbscan results/suggestive_strs_outliers.tsv \
  --out results/covid_suggestive_genes_with_outlier_STRs.tsv

echo "=== 5/5 cross-validation com LitCovid (se existir) ==="
if [ -f ../litcovid_validation/results/gene_literature_summary.tsv ]; then
  python3 crossvalidate.py \
    --litcovid-summary ../litcovid_validation/results/gene_literature_summary.tsv \
    --subset results/covid_suggestive_genes_with_outlier_STRs.tsv \
    --out results/crossvalidated_genes.tsv
else
  echo "Aviso: resultados do LitCovid ausentes; rode litcovid_validation antes para gerar crossvalidated_genes.tsv"
fi
echo "Pronto."
