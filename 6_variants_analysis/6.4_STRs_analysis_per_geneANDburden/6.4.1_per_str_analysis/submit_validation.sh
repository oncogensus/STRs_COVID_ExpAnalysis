#!/bin/bash
# submit_validation.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete a validacao cruzada DBSCAN<->LitCovid em cadeia via PBS, com
# dependencia afterok: a Parte 2 (dbscan_subset) so inicia apos a Parte 1
# (litcovid_validation) terminar com sucesso, garantindo que o
# gene_literature_summary.tsv exista para o crossvalidate.
#
#   bash submit_validation.sh
#
# Para submeter manualmente, um a um:
#   cd litcovid_validation && qsub litcovid.pbs
#   cd ../6.4.1.2_dbscan_subset && qsub -W depend=afterok:<JOB1> 2_run_dbscan_subset.pbs
cd "$(dirname "$0")"

cd litcovid_validation
L=$(qsub litcovid.pbs)
echo "litcovid_validation -> $L"

cd ../6.4.1.2_dbscan_subset
D=$(qsub -W depend=afterok:$L 2_run_dbscan_subset.pbs)
echo "dbscan_subset       -> $D"

echo
echo "Acompanhe com:  qstat -u $USER"
echo "Resultados: litcovid_validation/results/ e 6.4.1.2_dbscan_subset/results/"
