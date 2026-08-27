#!/bin/bash
# submit_pathway_crossvalidation.sh
# Submete o cross-validation funcional de vias (KEGG+Reactome) no env r_enrich_env.
cd "$(dirname "$0")"
qsub pathway_crossvalidation.pbs
