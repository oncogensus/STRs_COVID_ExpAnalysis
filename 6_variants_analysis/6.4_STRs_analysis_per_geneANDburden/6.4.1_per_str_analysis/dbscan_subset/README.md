# dbscan_subset  (Parte 2 da validação DBSCAN↔LitCovid)

Dados os **genes sugestivos do COVID-19 HG r7** (p < 1e-5, fenótipos A2/B2/C2,
`leave_23andme`), localizamos os STR loci da coorte que caem nesses genes e
**re-rodamos o DBSCAN** (mesmos parâmetros de `5_dbscan/outliers_search/5.2_dbscan_str.r`)
sobre esse subset, buscando outliers. Responde a: "nos genes sugestivos do COVID-19 HG,
a coorte tem STRs que o DBSCAN marca como outliers?"

## Etapas
1. `run_dbscan_subset.sh` orquestra:
   - `../extract_covid_genes.py` (reutilizado) → `results/covid_genes_suggestive.tsv` (p<1e-5).
   - `overlap_suggestive_strs.py` → `results/suggestive_gene_strs.tsv` + `data/suggestive_strs_ids.txt`.
   - subset de `5_dbscan/norm_test/STRs_normalized_residuals.tsv` (por STRs_ID) → `run_dbscan_subset.R` → `results/suggestive_strs_outliers.tsv`.
   - `summarize_subset.py` (filtra n_outliers>0) → `results/covid_suggestive_genes_with_outlier_STRs.tsv`.
   - `crossvalidate.py` (se existir o resumo do LitCovid) → `results/crossvalidated_genes.tsv`.

## Como rodar (cluster, env `igv`)
```bash
bash run_dbscan_subset.sh          # ou: qsub dbscan_subset.pbs  (submit_dbscan_subset.sh)
```
Para submeter as DUAS etapas (Parte 1 -> Parte 2) encadeadas via PBS, use na pasta
pai: `bash submit_validation.sh` (a Parte 2 depende da Parte 1 para o crossvalidate).
`CATALOG` e `RESIDUALS` têm defaults relativos à raiz do repo (3 níveis acima); sobrescreva se precisar.

## Pré-requisitos
- `genes.hg38.bed` e `data/*.leave_23andme_20220403.tsv.gz` já existem em
  `../` (rode `download_data.sh` + `build_gene_bed.py` em `covid19hg_evaluation/` antes, se não).
- `5_dbscan/norm_test/STRs_normalized_residuals.tsv` (saída de `5_dbscan`).

## Saídas
- `results/covid_genes_suggestive.tsv` — genes sugestivos (p<1e-5).
- `results/suggestive_gene_strs.tsv` — gene × STR da coorte nesse gene.
- `results/suggestive_strs_outliers.tsv` — outliers DBSCAN por STR (subset).
- `results/covid_suggestive_genes_with_outlier_STRs.tsv` — pares gene sugestivo × STR outlier.
- `results/crossvalidated_genes.tsv` — interseção com LitCovid (se disponível).
