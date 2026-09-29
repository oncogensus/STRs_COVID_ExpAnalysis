# 6.5 — Ancestry Analysis

Ancestry analysis: correlation between ancestry proportions (EthSEQ) and DBSCAN outlier metrics.

## Structure

```
6.5_ancestry_analysis/
├── 1_ancestry_comparison_cat.r         # Categorical comparison (Kruskal-Wallis + Dunn)
├── 2_ancestry_comparison_high_resolution.r  # High-resolution correlation (Spearman)
└── 3_ancestry_dataviz.ipynb             # Publication-ready visualization
```

## Data Flow

```
samples/STRs_analysis_dataset.tsv (unified dataset, stage 6.1)
    ↓
1_ancestry_comparison_cat.r  →  results/categorical_data/
    ↓                           (Kruskal-Wallis + Dunn post-hoc tests)
    ↓
2_ancestry_comparison_high_resolution.r  →  results/high_resolution/
    ↓                                       (continuous Spearman correlation)
    ↓
3_ancestry_dataviz.ipynb  →  results/high_resolution/dataviz/
                             (GT tables, heatmaps, ridge plots)
```

---

## 1. Categorical Comparison (`1_ancestry_comparison_cat.r`)

Compares allele distributions and DBSCAN outlier burden between categorical populations via Kruskal-Wallis and Dunn post-hoc tests.

**Input**: `STRs_analysis_dataset.tsv`

**QC Filters** (7,141 STRs passing filters):
- `n_clusters > 0`
- `noise_ratio <= 0.10`
- `n_outliers >= 1`

**Outputs** (`results/categorical_data/`):

| File | Description |
|---|---|
| `alleles_distribution_summary.csv` | Allele distribution summary |
| `alleles_kruskal_results.csv` | Kruskal-Wallis results (alleles) |
| `alleles_dunn_results.csv` | Dunn post-hoc results (alleles) |
| `dbscan_distribution_summary.csv` | DBSCAN distribution summary |
| `dbscan_kruskal_results.csv` | Kruskal-Wallis results (DBSCAN) |
| `dbscan_dunn_results.csv` | Dunn post-hoc results (DBSCAN) |
| `plotdata_alleles_long.csv` | Plot data (alleles, long) |
| `plotdata_alleles_wide.csv` | Plot data (alleles, wide) |
| `plotdata_dbscan_long.csv` | Plot data (DBSCAN, long) |
| `plotdata_dbscan_wide.csv` | Plot data (DBSCAN, wide) |
| `dbscan_qc_flags.csv` | DBSCAN QC flags |

**Execution**:
```bash
cd 6.5_ancestry_analysis
Rscript 1_ancestry_comparison_cat.r
```

---

## 2. High-Resolution Correlation (`2_ancestry_comparison_high_resolution.r`)

Correlates continuous ancestry proportions (EthSEQ) with DBSCAN outlier metrics (proportion and strength) per genomic region using Spearman correlation.

**Input**: `STRs_analysis_dataset.tsv`

**Outputs** (`results/high_resolution/`):

| File | Description |
|---|---|
| `plotdata_region_sample.csv` | Outlier metrics per region x sample |
| `correlation_full.csv` | Spearman rho and adjusted p-value |
| `ancestry_region_distribution_wide.csv` | Ancestry proportions per region |

**Execution**:
```bash
Rscript 2_ancestry_comparison_high_resolution.r
```

---

## 3. Visualization (`3_ancestry_dataviz.ipynb`)

Generates publication-ready tables and heatmaps of the ancestry results.

**Inputs**: CSVs from steps 1 and 2

**Outputs** (`results/high_resolution/dataviz/`):

| File | Description |
|---|---|
| `genomic_summary_per_region.html` | Regional table (ancestry + outlier metrics) |
| `heatmap_correlation_outlier_prop.png` | Spearman rho heatmap (outlier proportion) |
| `heatmap_correlation_outlier_strength.png` | Spearman rho heatmap (outlier strength) |
| `comprehensive_correlation_table.html` | Full table with significance |
| Ridge plots, boxplots | Allele distributions per ancestry |

**Execution**:
```bash
# Open in Jupyter
jupyter notebook 3_ancestry_dataviz.ipynb
```

---

## Environment

- `r_enrich_env` (micromamba): `data.table`, `rstatix`, `dplyr`, `gt`, `ggplot2` (pinned as `r_enrich_env.yaml`, repo root)