# 6.2 — Descriptive Analysis & Genome Visualization

Descriptive analysis and genomic visualization of the cohort STRs.

## Structure

```
6.2_desc_data_viz/
├── dataviz/
│   └── 1_genome_viz.ipynb        # Genomic visualization (Karyotype, Ideogram, Circos)
│
└── desc_analysis/
    ├── 1_desc_analysis_submit.r  # Submits all R scripts via PBS
    ├── 2_desc_analysis.ipynb     # Interactive notebook (exploratory analysis)
    ├── 3_str_coverage_per_patient.R   # STR coverage per patient
    ├── 4_str_coverage_gt_global.R     # GT table of global coverage
    ├── 5_dbscan_validation.R          # DBSCAN technical validation
    ├── 6_outlier_report.R             # Unified outlier report
    ├── 7_merged_table_allele_stats.R  # Allele statistics table
    ├── genome.txt                # Chromosome sizes (chr1-22, X, Y)
    └── desc_analysis_strs.log    # Execution log
```

---

## dataviz — Genome Visualization

### `1_genome_viz.ipynb`

Generates genomic visualizations using `regioneR` and `ggbio`:
- **Karyotype plot**: chromosomal distribution of the STRs
- **Ideogram**: chromosome bands with overlaid STRs
- **Circos plot**: multi-chromosome circular view

**Input**: `STRs_analysis_dataset.tsv` (unified dataset)
**Environment**: `r_viz` (micromamba)

---

## desc_analysis — Descriptive Analysis

### Genomic Coverage

#### `3_str_coverage_per_patient.R`

Computes, per patient, the STR genomic coverage after QC:
- `n_strs`: number of distinct STR loci
- `bp_covered`: covered bases (merged overlapping intervals per chr)
- `genome_frac`: fraction of the genome covered

**Inputs**: `--dataset STRs_analysis_dataset.tsv`, `--genome genome.txt`
**Outputs**: `coverage_by_patient.tsv`, `coverage_summary.tsv`

#### `4_str_coverage_gt_global.R`

Publication-ready GT table with the global coverage summary (reads `coverage_summary.tsv` from the previous step).

**Input**: `--summary coverage_summary.tsv`
**Outputs**: `table_coverage_global.html`, `table_coverage_global.csv`

---

### DBSCAN Validation

#### `5_dbscan_validation.R`

Dual-panel technical validation of the DBSCAN (all cohort loci):
- **Panel A**: Genotype distribution (1 Cluster, 2 Clusters, 3+, Unknown)
- **Panel B**: Noise tiers (High Quality < 0.10, Acceptable < 0.25, Other)

**Input**: `--str-catalog STRs_analysis_dataset.tsv`
**Outputs**: `dbscan_dual_panel_validation.png`, `quality_funnel_summary.csv`

#### `6_outlier_report.R`

Unified DBSCAN outlier report for the global cohort:
- General Summary (total, signal, no_signal)
- Cluster Distribution (1, 2, 3+ clusters)
- Noise Tiers (0-5%, 5-10%, 10-20%, 20-50%, >50%)
- Outlier Frequency (0 vs >0 outliers)

**Input**: `--str-catalog STRs_analysis_dataset.tsv`
**Outputs**: `unified_binary_outlier_report.csv`, `unified_technical_report.csv`, `Technical_Validation_Table.html`

---

### Allele Statistics

#### `7_merged_table_allele_stats.R`

Transposed (wide) table with descriptive statistics of `allele2_est` per genomic region + overall, comparing Case vs Control.

**Input**: `results/strs_by_locus_combo.csv`
**Outputs**: `table_allele_stats_merged.html`, `table_allele_stats_merged.csv`

---

## Execution

```bash
# Submit all R scripts via PBS
cd 6.2_desc_data_viz/desc_analysis
Rscript 1_desc_analysis_submit.r

# Or individually
Rscript 3_str_coverage_per_patient.R --dataset ../../samples/STRs_analysis_dataset.tsv --out-dir results/coverage
Rscript 4_str_coverage_gt_global.R --summary results/coverage/coverage_summary.tsv --out-dir results/coverage
Rscript 5_dbscan_validation.R --str-catalog ../../samples/STRs_analysis_dataset.tsv --out-dir results/dbscan
Rscript 6_outlier_report.R --str-catalog ../../samples/STRs_analysis_dataset.tsv --out-dir results/outlier
Rscript 7_merged_table_allele_stats.R
```

## Environment

- `r_enrich_env` (micromamba): `data.table`, `dplyr`, `gt`, `ggplot2`, `patchwork` (pinned as `r_enrich_env.yaml`, repo root)
- `r_viz` (micromamba): `regioneR`, `ggbio` (for dataviz) (pinned as `r_viz.yaml`, repo root)