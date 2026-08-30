# 6.4.1 — Per-STR Analysis

Análises dirigidas a genes/STRs específicos, combinando evidências de
**associação genética (COVID-19 HG r7)**, **outliers DBSCAN** e **vias
biológicas**.

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
3. `3_summarize_subset.py` — filtra n_outliers > 0 e anexa anotação do catálogo
   (`gene_id`, `gene_name`, `region`, `type`, `pops`) → `results/covid_suggestive_genes_with_outlier_STRs.tsv`
4. `4_annotate_catalog.pbs` — `annotate_catalog.py` enriquece resíduos (modo `sample`) e
   outliers (modo `str`) com o catálogo completo da coorte.

**Anotação do catálogo** (`annotate_catalog.py`, em `4_annotate_catalog.pbs`):
junta `samples/STRs_analysis_dataset.tsv` (via `$CATALOG`) trazendo **TODAS** as
colunas do catálogo **exceto** as métricas DBSCAN globais
(`n_outliers`, `outlier_samples`, `outlier_residuals`, `n_clusters`, `noise_ratio`).
- Modo `sample` (resíduos, tem `sample_id`): acrescenta `group`, `age`, `sex`,
  `allele1_est`, `allele2_est`, `depth`, `repeat_unit`, `gene_id`, `gene_name`,
  `region`, `chrom`, `start`, `end`, `pop` (ancestralidade), `contribution`, `type`.
- Modo `str` (outliers, só `STRs_ID`): acrescenta só colunas de locus
  (`repeat_unit`, `gene_id`, `gene_name`, `region`, `chrom`, `start`, `end`, `type`),
  preservando as métricas DBSCAN do subset.

**Arquivos gerados (enriquecidos com resíduos, métricas DBSCAN e catálogo da coorte)**:
- `data/suggestive_strs_residuals.tsv` — subset de `STRs_normalized_residuals.tsv`
  (um residuo `allele2_residuals` por amostra/STR) + anotação completa do catálogo
  por amostra (ancestralidade `pop`, `gene_id`, `gene_name`, `region`, `type`, etc.).
- `results/suggestive_strs_outliers.tsv` — 1 linha por STR com métricas DBSCAN
  (`n_samples`, `n_samples_valid`, `n_outliers`, `outlier_samples`,
  `outlier_residuals`, `n_clusters`, `noise_ratio`, `eps`, `minPts`, `cutoff`,
  `max_residual`, `mean_residual`) + anotação de locus do catálogo.
- `results/covid_suggestive_genes_with_outlier_STRs.tsv` — gene sugestivo × STR
  outlier: localização completa (gene + STR), todas as métricas/resíduos acima e
  anotação do catálogo (`gene_id`, `gene_name`, `region`, `type`, `pops`).

**Submissão PBS**: `submit_dbscan_subset.sh`.

---

### 6.4.1.3 — Cross-validation de vias (`6.4.1.3_pathway_crossvalidation/`)

Enriquecimento de vias biológicas nos genes outlier, com validação cruzada
por sub-amostragem.

**Ordem de execução**:
1. `1_pathway_crossvalidation.R` — enriquecimento de vias
2. `2_gene_crossvalidation.R` — validação cruzada por gene
3. `3_plot_pathways.R` — visualização dos resultados

**Submissão PBS**: `submit_pathway_crossvalidation.sh`, `submit_gene_crossvalidation.sh`.

---

## Ambiente

- `igv` (micromamba): scripts Python e download.
- `dbscan-r` (micromamba): `2_run_dbscan_subset.R`.
