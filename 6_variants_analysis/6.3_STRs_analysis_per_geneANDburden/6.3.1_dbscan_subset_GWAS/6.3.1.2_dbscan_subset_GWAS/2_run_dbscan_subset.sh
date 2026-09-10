#!/usr/bin/env bash
# run_dbscan_subset.sh  (Parte 2)
# Genes sugestivos (p<1e-5) e significativos (p<5e-8) COVID-19 HG -> STRs na coorte -> re-roda DBSCAN no subset.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
ROOT="$(cd "$HERE/../../.." && pwd)"
mkdir -p data results

# Caminhos do micromamba (a funcao `micromamba` do hook nao e herdada pelo subshell)
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/storage2/matheusbomfim/projects/micromamba}"
MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
DBSCAN_RSCRIPT="$MAMBA_ROOT_PREFIX/envs/dbscan-r/bin/Rscript"

# Garante data.table/dbscan no env dbscan-r (no-op se ja existirem)
"$MAMBA_BIN" install -n dbscan-r -c conda-forge -y r-data.table r-dbscan 2>&1 | tail -5 || \
  echo "AVISO: nao foi possivel instalar via micromamba; se os pacotes ja existirem no env dbscan-r, pode ignorar."
CATALOG="${CATALOG:-$ROOT/samples/STRs_analysis_dataset.tsv}"
RESIDUALS="${RESIDUALS:-$ROOT/5_global_dbscan/norm_test/STRs_normalized_residuals.tsv}"

echo "=== 1a/4 extrai genes sugestivos COVID-19 HG (p<1e-5) ==="
python3 ../extract_covid_genes.py \
  --sumstats \
    ../data/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz \
  --gene-bed ../genes.hg38.bed \
  --out results/covid_genes_suggestive.tsv \
  --p-thresh 1e-5

echo "=== 1b/4 extrai genes significativos COVID-19 HG (p<5e-8) ==="
python3 ../extract_covid_genes.py \
  --sumstats \
    ../data/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz \
    ../data/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz \
  --gene-bed ../genes.hg38.bed \
  --out results/covid_genes_significant.tsv \
  --p-thresh 5e-8

echo "=== 2/4 overlap genes sugestivos x STRs da coorte (por gene_name) ==="
python3 overlap_suggestive_strs.py \
  --catalog "$CATALOG" \
  --covid-genes results/covid_genes_suggestive.tsv \
  --out results/STRs_analysis_dataset_with_GWAS.tsv \
  --out-overlap results/suggestive_gene_strs.tsv \
  --ids-out data/suggestive_strs_ids.txt

echo "=== 3/4 subset de residuos + DBSCAN ==="
awk -F'\t' '
NR==FNR { ids[$1]; next }
FNR==1  { for(i=1;i<=NF;i++) if($i=="STRs_ID") c=i; print; next }
$(c) in ids
' data/suggestive_strs_ids.txt "$RESIDUALS" > data/suggestive_strs_residuals.tsv
# R roda no env micromamba dbscan-r (tem data.table/dbscan)
"$DBSCAN_RSCRIPT" run_dbscan_subset.R data/suggestive_strs_residuals.tsv results/suggestive_strs_outliers.tsv

echo "=== 4/4 junta outliers + resumo comparativo ==="
python3 3_summarize_subset.py \
  --unified results/STRs_analysis_dataset_with_GWAS.tsv \
  --dbscan results/suggestive_strs_outliers.tsv \
  --covid-genes-sig results/covid_genes_significant.tsv \
  --out results/covid_suggestive_genes_with_outlier_STRs.tsv \
  --summary results/summary_overlap_dbscan.tsv

echo "Pronto."
