#!/bin/bash
# submit_dbscan_subset.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete as 4 etapas de 6.4.1.2 em cadeia via PBS, com dependencia afterok:
#   1_overlap -> 2_dbscan -> 3_summarize -> 4_crossvalidate
#
#   bash submit_dbscan_subset.sh
cd "$(dirname "$0")"

D=$(qsub 1_overlap_suggestive_strs.pbs)
echo "1_overlap_suggestive_strs -> $D"

B=$(qsub -W depend=afterok:$D 2_run_dbscan_subset.pbs)
echo "2_run_dbscan_subset       -> $B"

S=$(qsub -W depend=afterok:$B 3_summarize_subset.pbs)
echo "3_summarize_subset        -> $S"

C=$(qsub -W depend=afterok:$S 4_crossvalidate.pbs)
echo "4_crossvalidate           -> $C"

echo
echo "Acompanhe com:  qstat -u $USER"
echo "Resultados em results/"
