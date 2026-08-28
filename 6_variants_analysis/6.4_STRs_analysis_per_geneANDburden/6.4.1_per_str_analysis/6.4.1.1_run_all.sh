#!/usr/bin/env bash
# run_all.sh
# Orquestra o pipeline COVID-19 HG x STRs da coorte.
#   bash run_all.sh
# Etapas:
#   1) download_data.sh   (baixa summary stats + GTF; pode pular se ja baixou)
#   2) build_gene_bed.py  (GTF -> genes.hg38.bed)
#   3) extract_covid_genes.py (summary -> covid_genes.tsv)
#   4) overlap_str_cohort.py  (covid_genes x catalogo STR -> covid_genes_with_cohort_STRs.tsv)
#
# Edite CATALOG para o caminho do catalogo completo de STRs da coorte.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# >>> AJUSTE: caminho do catalogo completo de STRs da coorte <<<
CATALOG="/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis/samples/STRs_analysis_dataset.tsv"
WINDOW=50000
P_THRESH=5e-8

echo "=== 1/4 download (pule se ja tiver ./data) ==="
bash download_data.sh || echo "download falhou/ignorado"

echo "=== 2/4 gene BED ==="
python3 build_gene_bed.py --gtf data/gencode.v39.primary_assembly.annotation.gtf.gz --out genes.hg38.bed

echo "=== 3/4 extrai genes COVID (p<$P_THRESH, janela $WINDOW) ==="
python3 extract_covid_genes.py \
  --sumstats \
     data/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz \
     data/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz \
     data/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz \
  --gene-bed genes.hg38.bed \
  --out covid_genes.tsv \
  --window "$WINDOW" --p-thresh "$P_THRESH"

echo "=== 4/4 cruza com STRs da coorte ==="
python3 overlap_str_cohort.py \
  --catalog "$CATALOG" \
  --covid-genes covid_genes.tsv \
  --out covid_genes_with_cohort_STRs.tsv \
  --window "$WINDOW"

echo
echo "Pronto. Veja covid_genes.tsv e covid_genes_with_cohort_STRs.tsv"
