#!/bin/bash
# submit_pathway_crossvalidation.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete as 4 etapas de 6.4.1.4 em cadeia via PBS:
#   0_outlier_genes -> 1_pathway -> 2_gene -> 3_plots
#
#   bash submit_pathway_crossvalidation.sh
cd "$(dirname "$0")"

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
