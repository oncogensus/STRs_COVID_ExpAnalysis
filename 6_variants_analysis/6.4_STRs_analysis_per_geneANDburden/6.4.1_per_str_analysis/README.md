# 6.4.1 — Per-STR Analysis

Análises dirigidas a genes/STRs específicos, combinando evidências de
**associação genética (COVID-19 HG r7)**, **outliers DBSCAN**, **literatura
COVID** e **vias biológicas**.

## Sub-etapas

### 6.4.1.1 — COVID-19 HG × STRs da coorte (`6.4.1.1_covid19hg_overlap/`)

Cruzamento entre genes significativos do COVID-19 Host Genetics Initiative (r7)
e o catálogo completo de STRs da coorte (STRling).

**Dados**: summary stats COVID-19 HG r7 (A2/B2/C2, `leave_23andme`), GENCODE v39,
catálogo de STRs (`STRs_analysis_dataset.tsv`).

**Método**: SNPs com `p < 5e-8` mapeados ao gene **somente dentro do corpo do
gene** (sem janela de flanco); cruzamento coord-a-coord com catálogo de STRs.

**Ordem de execução**:
1. `1_download_data.sh` — baixa summary stats + GTF
2. `2_build_gene_bed.py` — GTF → `genes.hg38.bed`
3. `3_extract_covid_genes.py` — summary → `covid_genes.tsv`
4. `4_overlap_str_cohort.py` — `covid_genes.tsv` × catálogo → `covid_genes_with_cohort_STRs.tsv`

**Orquestração**: `run_all.sh` (sequencial) ou `submit_all.sh` (PBS encadeado).

---

### 6.4.1.2 — DBSCAN subset — cross-validation (`6.4.1.2_dbscan_subset/`)

Re-execução do DBSCAN (mesmos parâmetros de `5_global_dbscan`) sobre STRs
localizados em genes sugestivos do COVID-19 HG r7 (`p < 1e-5`).

**Pré-requisito**: `6.4.1.1_covid19hg_overlap/` já rodado (precisa de
`genes.hg38.bed` e `data/*.tsv.gz`); `5_global_dbscan/norm_test/STRs_normalized_residuals.tsv`.

**Ordem de execução**:
1. `1_overlap_suggestive_strs.py` — gene × STR sugestivo → `results/suggestive_gene_strs.tsv`
2. `2_run_dbscan_subset.sh` / `2_run_dbscan_subset.R` — DBSCAN no subset → `results/suggestive_strs_outliers.tsv`
3. `3_summarize_subset.py` — filtra n_outliers > 0 → `results/covid_suggestive_genes_with_outlier_STRs.tsv`
4. `4_crossvalidate.py` — interseção com LitCovid → `results/crossvalidated_genes.tsv`

**Arquivos gerados (enriquecidos com resíduos, métricas DBSCAN e localização genômica)**:
- `data/suggestive_strs_residuals.tsv` — subset de `STRs_normalized_residuals.tsv`
  (um residuo `allele2_residuals` por amostra/STR) + colunas de localização
  (`str_chrom`, `str_start`, `str_end`, `repeatunit`, `str_locus`).
- `results/suggestive_strs_outliers.tsv` — 1 linha por STR com métricas DBSCAN
  (`n_samples`, `n_samples_valid`, `n_outliers`, `outlier_samples`,
  `outlier_residuals`, `n_clusters`, `noise_ratio`, `eps`, `minPts`, `cutoff`,
  `max_residual`, `mean_residual`) + localização genômica do STR.
- `results/covid_suggestive_genes_with_outlier_STRs.tsv` — gene sugestivo × STR
  outlier: traz localização completa (gene + STR) e todas as métricas/resíduos
  acima. `annotate_location.py` adiciona a localização aos arquivos de residuos/outliers.

**Submissão PBS**: `submit_dbscan_subset.sh`.

---

### 6.4.1.3 — Validação por literatura — LitCovid (`6.4.1.3_litcovid_validation/`)

Extração de literatura COVID que menciona genes dos STRs marcados como
outliers pelo DBSCAN global.

**Pré-requisito**: `5_global_dbscan/outliers_search/results_dbscan/outliers_per_str.tsv`.

**Ordem de execução**:
1. `1_get_outlier_genes.py` — catálogo → `data/outlier_genes.txt`
2. `2_download_litcovid.sh` — baixa `litcovid2pubtator.json.gz` (~2.3 GB)
3. `3_gene_literature.py` — mapeia gene → artigos → `results/gene_literature_summary.tsv`

**Submissão PBS**: `submit_litcovid.sh` ou `run_litcovid.sh`.

---

### 6.4.1.4 — Cross-validation de vias (`6.4.1.4_pathway_crossvalidation/`)

Enriquecimento de vias biológicas nos genes outlier, com validação cruzada
por sub-amostragem.

**Ordem de execução**:
1. `1_pathway_crossvalidation.R` — enriquecimento de vias
2. `2_gene_crossvalidation.R` — validação cruzada por gene
3. `3_plot_pathways.R` — visualização dos resultados

**Submissão PBS**: `submit_pathway_crossvalidation.sh`, `submit_gene_crossvalidation.sh`.

---

## Orquestrador de validação cruzada

`submit_validation.sh` (na raiz) encadeia as etapas 6.4.1.3 → 6.4.1.2 via PBS com
dependência `afterok`:

```bash
bash submit_validation.sh
# litcovid_validation (Parte 1) → dbscan_subset (Parte 2)
```

## Ambiente

- `igv` (micromamba): scripts Python e download.
- `dbscan-r` (micromamba): `2_run_dbscan_subset.R`.
