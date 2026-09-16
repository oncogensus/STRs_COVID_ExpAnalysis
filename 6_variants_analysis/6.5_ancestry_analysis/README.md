# 6.5 — Ancestry Analysis

Análise de ancestralidade: correlação entre proporções de ancestralidade (EthSEQ) e métricas DBSCAN outliers.

## Estrutura

```
6.5_ancestry_analysis/
├── 1_ancestry_comparation_cat.r         # Comparação categórica (Kruskal-Wallis + Dunn)
├── 2_ancestry_comparation_high_resolution.r  # Correlação alta resolução (Spearman)
└── 3_ancestry_dataviz.ipynb             # Visualização publication-ready
```

## Fluxo de Dados

```
samples/STRs_analysis_dataset.tsv (dataset unificado, etapa 6.1)
    ↓
1_ancestry_comparation_cat.r  →  results/categorical_data/
    ↓                           (testes Kruskal-Wallis + Dunn post-hoc)
    ↓
2_ancestry_comparation_high_resolution.r  →  results/high_resolution/
    ↓                                       (correlação Spearman contínua)
    ↓
3_ancestry_dataviz.ipynb  →  results/high_resolution/dataviz/
                             (tabelas GT, heatmaps, ridge plots)
```

---

## 1. Comparação Categórica (`1_ancestry_comparation_cat.r`)

Compara distribuições de alleles e burden de outliers DBSCAN entre populações categóricas via Kruskal-Wallis e Dunn post-hoc.

**Input**: `STRs_analysis_dataset.tsv`

**Filtros QC**:
- `n_clusters > 0`
- `noise_ratio <= 0.10`
- `n_outliers >= 1`

**Outputs** (`results/categorical_data/`):

| Arquivo | Descrição |
|---|---|
| `alleles_distribution_summary.csv` | Resumo de distribuição de alleles |
| `alleles_kruskal_results.csv` | Resultados Kruskal-Wallis (alleles) |
| `alleles_dunn_results.csv` | Resultados Dunn post-hoc (alleles) |
| `dbscan_distribution_summary.csv` | Resumo de distribuição DBSCAN |
| `dbscan_kruskal_results.csv` | Resultados Kruskal-Wallis (DBSCAN) |
| `dbscan_dunn_results.csv` | Resultados Dunn post-hoc (DBSCAN) |
| `plotdata_alleles_long.csv` | Dados para plotagem (alleles, long) |
| `plotdata_alleles_wide.csv` | Dados para plotagem (alleles, wide) |
| `plotdata_dbscan_long.csv` | Dados para plotagem (DBSCAN, long) |
| `plotdata_dbscan_wide.csv` | Dados para plotagem (DBSCAN, wide) |
| `dbscan_qc_flags.csv` | Flags de QC do DBSCAN |

**Execução**:
```bash
cd 6.5_ancestry_analysis
Rscript 1_ancestry_comparation_cat.r
```

---

## 2. Correlação Alta Resolução (`2_ancestry_comparation_high_resolution.r`)

Correlaciona proporções contínuas de ancestralidade (EthSEQ) com métricas de outliers DBSCAN (proporção e força) por região genômica usando Spearman correlation.

**Input**: `STRs_analysis_dataset.tsv`

**Outputs** (`results/high_resolution/`):

| Arquivo | Descrição |
|---|---|
| `plotdata_region_sample.csv` | Métricas de outlier por região × sample |
| `correlation_full.csv` | Spearman rho e p-valor ajustado |
| `ancestry_region_distribution_wide.csv` | Proporções de ancestralidade por região |

**Execução**:
```bash
Rscript 2_ancestry_comparation_high_resolution.r
```

---

## 3. Visualização (`3_ancestry_dataviz.ipynb`)

Gera tabelas publication-ready e heatmaps dos resultados de ancestralidade.

**Inputs**: CSVs das etapas 1 e 2

**Outputs** (`results/high_resolution/dataviz/`):

| Arquivo | Descrição |
|---|---|
| `genomic_summary_per_region.html` | Tabela regional (ancestry + outlier metrics) |
| `heatmap_correlation_outlier_prop.png` | Heatmap Spearman rho (proporção outlier) |
| `heatmap_correlation_outlier_strength.png` | Heatmap Spearman rho (força outlier) |
| `comprehensive_correlation_table.html` | Tabela completa com significância |
| Ridge plots, boxplots | Distribuições de allele por ancestralidade |

**Execução**:
```bash
# Abrir no Jupyter
jupyter notebook 3_ancestry_dataviz.ipynb
```

---

## Ambiente

- `r_enrich_env` (micromamba): `data.table`, `rstatix`, `dplyr`, `gt`, `ggplot2`
