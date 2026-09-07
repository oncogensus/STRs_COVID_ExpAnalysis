#!/bin/bash
# run_all.sh -- roda localmente (login node) o modulo 6.4.7 em sequencia.
# Os passos 1 e 2 rodam para gwas_burden e rna_burden; o 3 usa os genes
# significativos; o 4 usa os arquivos de DEGs (edite DEG_DIR no .pbs).
set -e
REPO="/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"
WORKDIR="$REPO/6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment"
MAMBA="/storage2/matheusbomfim/projects/micromamba"

export MAMBA_ROOT_PREFIX="$MAMBA"
eval "$(/storage2/matheusbomfim/projects/micromamba/bin/micromamba shell hook --shell bash)"
micromamba activate r_enrich_env
R_BIN="$(command -v Rscript)"
[ -x "$R_BIN" ] || { echo "ERRO: Rscript ausente"; exit 1; }

cd "$WORKDIR"
git stash; git pull; git stash pop || true

for STRAT in gwas_burden rna_burden; do
  echo "===== 1. SKAT-O por gene ($STRAT) ====="
  "$R_BIN" 1_skat_o_per_gene.R --strategy "$STRAT"
done

for STRAT in gwas_burden rna_burden; do
  echo "===== 2. SKAT-O por via ($STRAT) ====="
  "$R_BIN" 2_skat_o_per_pathway.R --strategy "$STRAT"
done

echo "===== 3. Enriquecimento de vias ====="
bash 3_run_enrichment_for_sources.sh || echo "AVISO: passo 3 parcial"

echo "===== 4. Enriquecimento DE x STR ====="
bash 4_de_enrichment.pbs 2>/dev/null || echo "AVISO: rode 4_de_enrichment.pbs via qsub (DEG_DIR) "

echo "=== FIM run_all (6.4.7) ==="