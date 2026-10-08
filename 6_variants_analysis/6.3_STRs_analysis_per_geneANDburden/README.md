# 6.3 — STRs Analysis per Gene & Burden

Main module for the analysis of COVID-19-associated STRs: cross-referencing against RNA-seq DEGs, visualization, statistical testing, and visual inspection via IGV.

## Sub-modules

| Sub-module | Description | Documentation |
|---|---|---|
| `6.3.1_pre_processing/` | Preprocessing and differential expression of public COVID-19 transcriptomic datasets (GEO) | [`6.3.1_pre_processing/README.md`](6.3.1_pre_processing/README.md) |
| `6.3.2_RNA_data_analysis/` | DEG x STR crossing, raincloud plots, Mann-Whitney U, no-overlap and outlier detail tables | [`6.3.2_RNA_data_analysis/README.md`](6.3.2_RNA_data_analysis/README.md) |
| `6.3.3_burden_analysis/` | Relative outlier burden: Mann-Whitney U + Firth logistic regression (global, DEGs, per region) | [`6.3.3_burden_analysis/README.md`](6.3.3_burden_analysis/README.md) |
| `6.3.4_igv_per_variant/` | IGV.js visual inspection of each STR with outliers | [`6.3.4_igv_per_variant/README.md`](6.3.4_igv_per_variant/README.md) |

## Data Flow

```
samples/STRs_analysis_dataset.tsv (unified dataset, stage 6.1)
    ↓
6.3.2: DEG x STR crossing (per GSE intervention)
    ↓
    intervention_strs.tsv / intervention_outliers.tsv / intervention_summary.tsv
    ↓
    raincloud + publication tables, Mann-Whitney U, no-overlap, outlier detail tables
    ↓
6.3.3: relative outlier burden (global, DEGs, per region) → burden_*.csv
    ↓
6.3.4: BEDs + BAM mapping → IGV.js per variant
```

## Prerequisites

- `6.1_merge_datasets/` must be run first (generates `STRs_analysis_dataset.tsv`)
- `5_global_dbscan/` must have been run (global DBSCAN outliers)
- GSE subdirectories with DEG TSVs (for 6.3.2)

## Execution Order

```bash
# 6.3.2: DEG x STR crossing and analyses
cd 6.3.2_RNA_data_analysis/6.3.2.2_RNA_matrix
qsub 1_cross_intervention_STRs.pbs
# wait for completion
qsub 2_plot_intervention_summary.pbs
qsub 3_submit_raincloud.pbs
qsub 4_mann_whitney_allele.pbs
qsub 5_no_overlap_analysis.pbs
qsub 6_outlier_detail_table.pbs  # note: PBS file is named 8_outlier_detail_table.pbs

# 6.3.3: burden analyses
cd ../6.3.3_burden_analysis
Rscript 1_burden_analysis.R

# 6.3.4: IGV visualization
cd ../../6.3.4_igv_per_variant
qsub 1_generate_beds.pbs
bash 3_run_all.sh
```

Detailed script-level documentation, inputs/outputs and environments are in the [6.3.2](6.3.2_RNA_data_analysis/README.md) and [6.3.4](6.3.4_igv_per_variant/README.md) sub-module READMEs.