#!/bin/bash
# submit_gene_crossvalidation.sh
# Submete o cross-validation de genes (agnostic x GWAS-filtered) no env r_enrich_env.
cd "$(dirname "$0")"
qsub gene_crossvalidation.pbs
