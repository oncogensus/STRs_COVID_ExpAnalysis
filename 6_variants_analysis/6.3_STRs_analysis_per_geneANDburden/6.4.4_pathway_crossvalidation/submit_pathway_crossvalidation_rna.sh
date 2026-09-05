#!/bin/bash
# submit_pathway_crossvalidation_rna.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete a cross-validation GWAS x RNA de 6.4.4 em cadeia via PBS:
#   1_pathway_rna -> 2_gene_rna -> 3_plots_rna
# Pre-requisito: cross_DEGs_STRs.py (6.4.2.2) gerou results/rna_outlier_genes.tsv
#
#   bash submit_pathway_crossvalidation_rna.sh
cd "$(dirname "$0")"

RNA_FILE="../6.4.2_dbscan_subset_RNA/6.4.2.2_RNA_matrix/results/rna_outlier_genes.tsv"
[ -f "$RNA_FILE" ] || { echo "ERRO: $RNA_FILE ausente (rode cross_DEGs_STRs.py em 6.4.2.2)"; exit 1; }

P=$(qsub 1_pathway_crossvalidation_rna.pbs)
echo "1_pathway_crossvalidation (GWAS x RNA) -> $P"

G=$(qsub -W depend=afterok:$P 2_gene_crossvalidation_rna.pbs)
echo "2_gene_crossvalidation    (GWAS x RNA)  -> $G"

V=$(qsub -W depend=afterok:$G 3_plot_pathways_rna.pbs)
echo "3_plot_pathways           (GWAS x RNA)  -> $V"

echo
echo "Acompanhe com:  qstat -u $USER"
echo "Resultados em results/ (prefixo gwas_rna_)"