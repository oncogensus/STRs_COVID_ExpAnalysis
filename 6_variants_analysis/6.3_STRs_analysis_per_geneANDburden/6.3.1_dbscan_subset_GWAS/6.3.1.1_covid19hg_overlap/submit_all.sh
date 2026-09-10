#!/bin/bash
# submit_all.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete o pipeline COVID-19 HG x STRs em cadeia, com dependencia
# afterok: cada etapa so inicia apos a anterior terminar com sucesso.
#
#   bash submit_all.sh
#
# Para submeter manualmente, um a um:
#   qsub download_data.pbs
#   qsub build_gene_bed.pbs
#   qsub extract_covid_genes.pbs
#   qsub overlap_str_cohort.pbs
cd "$(dirname "$0")"

D=$(qsub 1_download_data.pbs)
echo "download        -> $D"
B=$(qsub -W depend=afterok:$D 2_build_gene_bed.pbs)
echo "build_gene_bed  -> $B"
E=$(qsub -W depend=afterok:$B 3_extract_covid_genes.pbs)
echo "extract_genes   -> $E"
O=$(qsub -W depend=afterok:$E 4_overlap_str_cohort.pbs)
echo "overlap         -> $O"
echo
echo "Acompanhe com:  qstat -u $USER"
echo "Resultados: covid_genes.tsv e covid_genes_with_cohort_STRs.tsv"
