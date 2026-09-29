# 6.3 — STRs Analysis per Gene & Burden

Main module for the analysis of COVID-19-associated STRs: cross-referencing against RNA-seq DEGs, visualization and inspection via IGV.

## Structure

```
6.3_STRs_analysis_per_geneANDburden/
├── .gitignore
├── 6.3.1_pre_processing/          # (reserved for future preprocessing)
├── 6.3.2_RNA_data_analysis/       # DEG x STR crossing + statistical analysis
│   ├── 6.3.2.2_RNA_matrix/
│   │   ├── 1_cross_intervention_STRs.py   # DEGs x STRs per intervention
│   │   ├── 2_plot_intervention_summary.R  # Raincloud + publication table
│   │   ├── 3_plot_raincloud_per_locus.R   # Raincloud per locus (outliers)
│   │   ├── 4_mann_whitney_allele.R        # Mann-Whitney U (case vs control)
│   │   ├── 5_no_overlap_analysis.R        # No-overlap analysis
│   │   └── 8_outlier_detail_table.R       # Detailed outlier table
│   └── README.md
│
└── 6.3.4_igv_per_variant/          # IGV visualization per variant
    ├── 1_generate_beds.R           # Generates BEDs + BAM mapping
    ├── 2_igv_variant.sh            # IGV.js for 1 STR
    ├── 3_run_all.sh                # Generates IGV scripts for all STRs
    ├── scripts/                    # Generated IGV scripts per variant
    └── README.md
```

## Data Flow

```
samples/STRs_analysis_dataset.tsv (unified dataset, stage 6.1)
    ↓
1_cross_intervention_STRs.py  →  intervention_strs.tsv (all STRs in DEG genes)
    ↓                            intervention_outliers.tsv (only STRs with DBSCAN outliers)
    ↓                            intervention_summary.tsv (summary per gene/intervention)
    ↓
2_plot_intervention_summary.R →  intervention_raincloud.png + publication tables
    ↓
3_plot_raincloud_per_locus.R  →  Raincloud per locus (DBSCAN outliers only)
    ↓
4_mann_whitney_allele.R       →  Mann-Whitney U (case vs control)
    ↓
5_no_overlap_analysis.R       →  Loci without case/control overlap
    ↓
8_outlier_detail_table.R      →  Detailed publication-ready table
    ↓
1_generate_beds.R             →  BEDs + BAM mapping
    ↓
2_igv_variant.sh              →  IGV.js for visual inspection
```

## Prerequisites

- `6.1_merge_datasets/` must be run first (generates `STRs_analysis_dataset.tsv`)
- `5_global_dbscan/` must have been run (global DBSCAN outliers)
- GSE subdirectories with DEG TSVs (for 6.3.2)

## Environment

| Script | Environment (micromamba) |
|---|---|
| `1_cross_intervention_STRs.py` | `str` |
| `2_plot_intervention_summary.R` | `r_enrich_env` |
| `3_plot_raincloud_per_locus.R` | `r_enrich_env` |
| `4_mann_whitney_allele.R` | `r_enrich_env` |
| `5_no_overlap_analysis.R` | `r_enrich_env` |
| `8_outlier_detail_table.R` | `r_enrich_env` |
| `1_generate_beds.R` | `igv` |

## Execution

```bash
# 6.3.2: DEG x STR crossing
cd 6.3.2_RNA_data_analysis/6.3.2.2_RNA_matrix
qsub 1_cross_intervention_STRs.pbs
# wait for completion
qsub 2_plot_intervention_summary.pbs
qsub 3_submit_raincloud.pbs
qsub 4_mann_whitney_allele.pbs
qsub 5_no_overlap_analysis.pbs
qsub 8_outlier_detail_table.pbs

# 6.3.4: IGV visualization
cd ../../6.3.4_igv_per_variant
qsub 1_generate_beds.pbs
bash 3_run_all.sh
```