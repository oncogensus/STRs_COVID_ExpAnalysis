#!/bin/bash
# submit_pathway_crossvalidation.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete as 4 etapas de 6.4.4 em cadeia via PBS:
#   0_outlier_genes -> 1_pathway -> 2_gene -> 3_plots
#
#   bash submit_pathway_crossvalidation.sh
cd "$(dirname "$0")"

GWAS_DIR="../6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS"

# --- Corrigir results/results/ duplicado (se existir) ---
if [ -d "$GWAS_DIR/results/results" ]; then
  echo "Corrigindo results/results/ duplicado..."
  mv "$GWAS_DIR/results/results/"* "$GWAS_DIR/results/" 2>/dev/null
  rmdir "$GWAS_DIR/results/results" 2>/dev/null
fi
if [ -d "$GWAS_DIR/data/data" ]; then
  echo "Corrigindo data/data/ duplicado..."
  mv "$GWAS_DIR/data/data/"* "$GWAS_DIR/data/" 2>/dev/null
  rmdir "$GWAS_DIR/data/data" 2>/dev/null
fi

# --- Pre-verificacao: resultados de 6.4.1.2 precisam existir ---
[ -f "$GWAS_DIR/results/covid_suggestive_genes_with_outlier_STRs.tsv" ] || \
  { echo "ERRO: rode 6.4.1.2_dbscan_subset_GWAS antes (results/covid_suggestive_genes_with_outlier_STRs.tsv ausente)"; exit 1; }

[ -f "$GWAS_DIR/data/suggestive_strs_residuals.tsv" ] || \
  { echo "ERRO: rode 6.4.1.2_dbscan_subset_GWAS antes (data/suggestive_strs_residuals.tsv ausente)"; exit 1; }

O=$(qsub 0_get_outlier_genes.pbs)
echo "0_get_outlier_genes      -> $O"

P=$(qsub -W depend=afterok:$O 1_pathway_crossvalidation.pbs)
echo "1_pathway_crossvalidation -> $P"

G=$(qsub -W depend=afterok:$P 2_gene_crossvalidation.pbs)
echo "2_gene_crossvalidation    -> $G"

V=$(qsub -W depend=afterok:$G 3_plot_pathways.pbs)
echo "3_plot_pathways           -> $V"

echo
echo "Acompanhe com:  qstat -u $USER"
echo "Resultados em results/"
