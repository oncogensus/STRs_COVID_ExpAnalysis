# 6.2 — Descriptive Analysis & Genome Visualization

Análise descritiva e visualização genômica dos STRs da coorte.

## Estrutura

```
6.2_desc_data_viz/
├── dataviz/
│   └── 1_genome_viz.ipynb        # Visualização genômica (Karyotype, Ideogram, Circos)
│
└── desc_analysis/
    ├── 1_desc_analysis_submit.r  # Submete todos os scripts R do PBS
    ├── 2_desc_analysis.ipynb     # Notebook interativo (análise exploratória)
    ├── 3_str_coverage_per_patient.R   # Cobertura STR por paciente
    ├── 4_str_coverage_gt_global.R     # Tabela GT da cobertura global
    ├── 5_dbscan_validation.R          # Validação técnica do DBSCAN
    ├── 6_outlier_report.R             # Relatório unificado de outliers
    ├── 7_merged_table_allele_stats.R  # Tabela de estatísticas de allele
    ├── genome.txt                # Tamanhos de cromossomo (chr1-22, X, Y)
    └── desc_analysis_strs.log    # Log de execução
```

---

## dataviz — Visualização Genômica

### `1_genome_viz.ipynb`

Gera visualizações genômicas usando `regioneR` e `ggbio`:
- **Karyotype plot**: distribuição cromossômica dos STRs
- **Ideogram**: bandas cromossômicas com STRs sobrepostos
- **Circos plot**: visão circular multi-cromossomo

**Input**: `STRs_analysis_dataset.tsv` (dataset unificado)
**Ambiente**: `r_viz` (micromamba)

---

## desc_analysis — Análise Descritiva

### Cobertura Genômica

#### `3_str_coverage_per_patient.R`

Calcula, por paciente, a cobertura genômica de STRs após QC:
- `n_strs`: número de loci STR distintos
- `bp_covered`: bases cobertas (intervalos sobrepostos mergeados por chr)
- `genome_frac`: fração do genoma coberta

**Inputs**: `--dataset STRs_analysis_dataset.tsv`, `--genome genome.txt`
**Outputs**: `coverage_by_patient.tsv`, `coverage_summary.tsv`

#### `4_str_coverage_gt_global.R`

Tabela GT publication-ready com resumo global da cobertura (lê `coverage_summary.tsv` do passo anterior).

**Input**: `--summary coverage_summary.tsv`
**Outputs**: `table_coverage_global.html`, `table_coverage_global.csv`

---

### Validação DBSCAN

#### `5_dbscan_validation.R`

Painel dual de validação técnica do DBSCAN (todos os loci da coorte):
- **Painel A**: Distribuição de genotipos (1 Cluster, 2 Clusters, 3+, Unknown)
- **Painel B**: Tiers de noise (High Quality < 0.10, Acceptable < 0.25, Other)

**Input**: `--str-catalog STRs_analysis_dataset.tsv`
**Outputs**: `dbscan_dual_panel_validation.png`, `quality_funnel_summary.csv`

#### `6_outlier_report.R`

Relatório unificado de outliers DBSCAN para a coorte global:
- General Summary (total, signal, no_signal)
- Cluster Distribution (1, 2, 3+ clusters)
- Noise Tiers (0-5%, 5-10%, 10-20%, 20-50%, >50%)
- Outlier Frequency (0 vs >0 outliers)

**Input**: `--str-catalog STRs_analysis_dataset.tsv`
**Outputs**: `unified_binary_outlier_report.csv`, `unified_technical_report.csv`, `Technical_Validation_Table.html`

---

### Estatísticas de Allele

#### `7_merged_table_allele_stats.R`

Tabela transposta (wide) com estatísticas descritivas de `allele2_est` por região genômica + overall, comparando Case vs Control.

**Input**: `results/strs_by_locus_combo.csv`
**Outputs**: `table_allele_stats_merged.html`, `table_allele_stats_merged.csv`

---

## Execução

```bash
# Submeter todos os scripts R via PBS
cd 6.2_desc_data_viz/desc_analysis
Rscript 1_desc_analysis_submit.r

# Ou individualmente
Rscript 3_str_coverage_per_patient.R --dataset ../../samples/STRs_analysis_dataset.tsv --out-dir results/coverage
Rscript 4_str_coverage_gt_global.R --summary results/coverage/coverage_summary.tsv --out-dir results/coverage
Rscript 5_dbscan_validation.R --str-catalog ../../samples/STRs_analysis_dataset.tsv --out-dir results/dbscan
Rscript 6_outlier_report.R --str-catalog ../../samples/STRs_analysis_dataset.tsv --out-dir results/outlier
Rscript 7_merged_table_allele_stats.R
```

## Ambiente

- `r_enrich_env` (micromamba): `data.table`, `dplyr`, `gt`, `ggplot2`, `patchwork`
- `r_viz` (micromamba): `regioneR`, `ggbio` (para dataviz)
