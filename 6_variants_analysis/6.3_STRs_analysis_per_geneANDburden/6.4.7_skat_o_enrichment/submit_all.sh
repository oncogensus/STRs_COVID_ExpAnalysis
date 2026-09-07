#!/bin/bash
# submit_all.sh -- submete os 4 passos do modulo 6.4.7 como jobs PBS.
# Uso: qsub submit_all.sh   (submete um job que encadeia os demais)
#PBS -q workq
#PBS -N skato_suite
#PBS -o skato_suite.out
#PBS -e skato_suite.err
#PBS -V

REPO="/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"
WORKDIR="$REPO/6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment"

cd "$WORKDIR" || exit 1

echo "==[1/4] SKAT-O por gene =="
qsub 1_skat_o_per_gene.pbs

echo "==[2/4] SKAT-O por via =="
qsub 2_skat_o_per_pathway.pbs

echo "==[3/4] Enriquecimento de vias (clusterProfiler) =="
qsub 3_enrichment_clusterprofiler.pbs

echo "==[4/4] Enriquecimento DE x STR =="
qsub 4_de_enrichment.pbs

echo "Jobs submetidos. Acompanhe com qstat."