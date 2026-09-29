# 6.1 — Merge Datasets

Merges genomic annotation, global DBSCAN results, ancestry and demographic data into a unified dataset.

## Objective

Consolidate the outputs from previous stages (GTF annotation, global DBSCAN, EthSEQ, phenotype) into a single `STRs_analysis_dataset.tsv` file that feeds all of stage 6.

## Structure

```
6.1_merge_datasets/
└── 1_merge_datasets.r      # Integration script
```

## Inputs

| File | Source | Description |
|---|---|---|
| `STRs_annotated_region.tsv` | `3_gtf_annot/samples/` | Genomic annotation per STR |
| `outliers_per_str.tsv` | `5_global_dbscan/outliers_search/results_dbscan/` | Global DBSCAN metrics per STR |
| `Report.txt` | `4_ancestry/EthSEQ_Results_3D/` | Ancestry assignment |
| `samples_infos.csv` | `samples/` | Demographic data (group, age, sex) |

## Output

`samples/STRs_analysis_dataset.tsv` — unified dataset with columns:

| Category | Columns |
|---|---|
| IDs & Clinical | `STRs_ID`, `group`, `age`, `sex` |
| STR Metrics | `allele1_est`, `allele2_est`, `depth` |
| Genomic Context | `repeat_unit`, `gene_id`, `gene_name`, `region`, `chrom`, `start`, `end`, `sample_id` |
| Population Genetics | `pop`, `contribution`, `type` |
| Global DBSCAN Outliers | `n_outliers_dbscan_global`, `outlier_samples_dbscan_global`, `outlier_residuals_dbscan_global`, `n_clusters_dbscan_global`, `noise_ratio_dbscan_global` |

## Execution

```bash
cd 6.1_merge_datasets
Rscript 1_merge_datasets.r
```

## Environment

- `r_enrich_env` (micromamba): `data.table`, `dplyr`, `stringr` (pinned as `r_enrich_env.yaml`, repo root)