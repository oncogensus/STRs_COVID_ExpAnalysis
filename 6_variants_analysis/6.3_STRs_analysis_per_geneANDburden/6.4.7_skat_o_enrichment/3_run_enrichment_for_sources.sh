#!/bin/bash
# 3_run_enrichment_for_sources.sh -- roda o enriquecimento de vias para
# múltiplas listas de genes-fonte (SKAT-O significativos, burden hits,
# outliers GWAS/RNA). Chamado pelo run_all.sh ou de um PBS.
set -e
REPO="/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"
WORKDIR="$REPO/6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment"
MAMBA="/storage2/matheusbomfim/projects/micromamba"

export MAMBA_ROOT_PREFIX="$MAMBA"
eval "$(/storage2/matheusbomfim/projects/micromamba/bin/micromamba shell hook --shell bash)"
micromamba activate r_enrich_env
R_BIN="$(command -v Rscript)"
[ -x "$R_BIN" ] || { echo "ERRO: Rscript ausente"; exit 1; }

cd "$WORKDIR"
mkdir -p results_enrichment

# extrai coluna 'gene' de um TSV e grava lista (arquivo, destino, label)
extract_gene_col() {
  local src="$1" dst="$2"
  if [ -f "$src" ]; then
    awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="gene")c=i; next} NR>1 && c{print $c}' \
      "$src" | sort -u > "$dst"
    echo "$(wc -l < "$dst") genes: $dst"
  else
    echo "AVISO: ausente: $src"
    : > "$dst"
  fi
}

run_enrich() {
  local genes="$1" label="$2"
  if [ ! -s "$genes" ]; then
    echo "-- Sem genes para $label --"
    return 0
  fi
  echo "===== Enriquecimento: $label ====="
  "$R_BIN" 3_enrichment_clusterprofiler.R --gene-file "$genes" \
    --label "$label" --out-dir "$WORKDIR/results_enrichment" || {
      echo "AVISO: enriquecimento $label falhou"; true; }
}

# Fontes: SKAT-O corrigido (BH)
extract_gene_col "$WORKDIR/results_skat_o_genes/skat_o_per_gene_hits_corrected.tsv" \
                 "$WORKDIR/results_enrichment/genes_skato_gwas_burden.txt"
extract_gene_col "$WORKDIR/results_skat_o_genes_rna/skat_o_per_gene_hits_corrected.tsv" \
                 "$WORKDIR/results_enrichment/genes_skato_rna_burden.txt"

# Fontes: burden hits (uncorrected) das duas estrategias
extract_gene_col "$WORKDIR/results_skat_o_genes/burden_per_gene_hits_uncorrected.tsv" \
                 "$WORKDIR/results_enrichment/genes_burden_gwas_unc.txt"
extract_gene_col "$WORKDIR/results_skat_o_genes_rna/burden_per_gene_hits_uncorrected.tsv" \
                 "$WORKDIR/results_enrichment/genes_burden_rna_unc.txt"

# Fontes: outlier genes (GWAS e RNA) -- gerados pelo comparativo 6.4.4
extract_gene_col \
  "$REPO/6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.4_pathway_crossvalidation/results_gwas_rna_comparison/outlier_genes_union.tsv" \
  "$WORKDIR/results_enrichment/genes_outlier_union.txt"

run_enrich "$WORKDIR/results_enrichment/genes_skato_gwas_burden.txt"  skato_gwas_burden
run_enrich "$WORKDIR/results_enrichment/genes_skato_rna_burden.txt"   skato_rna_burden
run_enrich "$WORKDIR/results_enrichment/genes_burden_gwas_unc.txt"    burden_gwas_unc
run_enrich "$WORKDIR/results_enrichment/genes_burden_rna_unc.txt"     burden_rna_unc
run_enrich "$WORKDIR/results_enrichment/genes_outlier_union.txt"      outlier_union

echo "=== FIM 3_run_enrichment_for_sources.sh ==="