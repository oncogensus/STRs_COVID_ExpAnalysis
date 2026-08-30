#!/bin/bash
# submit_gene_crossvalidation.sh  (roda NO NODE DE LOGIN, nao como PBS job)
# Submete somente a etapa 2 (gene cross-validation).
# Para rodar tudo em cadeia, use submit_pathway_crossvalidation.sh.
#
#   bash submit_gene_crossvalidation.sh
cd "$(dirname "$0")"
qsub 2_gene_crossvalidation.pbs
