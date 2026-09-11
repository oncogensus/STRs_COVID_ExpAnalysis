# 6.3.2 — RNA-seq Analysis

Cross-referencing differentially expressed genes (DEGs) from public RNA-seq datasets (GSE157103, GSE188847, GSE183533) with the cohort STR catalog, with DBSCAN outlier annotation.

## Structure

```
6.3.2_RNA_data_analysis/
├── 6.3.2.1_descriptive_analysis/    Descriptive analysis of RNA-seq outliers
│   ├── compare_gwas_rna.R           STR x patient overlap, case vs control stats
│   └── compare_gwas_rna.pbs
│
└── 6.3.2.2_RNA_matrix/              Per-intervention DEG x STR crossing + visualization
    ├── 1_cross_intervention_STRs.py   DEGs x STRs crossing (per intervention)
    ├── 1_cross_intervention_STRs.pbs
    ├── 2_plot_intervention_summary.R  Raincloud + publication table (per intervention)
    ├── 2_plot_intervention_summary.pbs
    ├── 3_plot_raincloud_per_locus.R   Raincloud per locus (outlier DBSCAN)
    └── 3_submit_raincloud.pbs
```

---

## 6.3.2.1 — Descriptive Analysis (`compare_gwas_rna.R`)

Descriptive analysis of RNA-seq outliers: STRs/genes with DBSCAN global outliers, patient-level allele distributions comparing case vs. control groups (no statistical tests).

**Inputs**:
- `intervention_strs.tsv` from `6.3.2.2_RNA_matrix/results/`
- `STRs_analysis_dataset.tsv` (cohort catalog)

**Outputs** (in `results/`):
| File | Description |
|---|---|
| `rna_outlier_sets.tsv` | STRs with outliers, gene, region, repeat_unit |
| `rna_genes_summary.tsv` | Genes: STR count, outlier count, GSE provenance |
| `patient_str.tsv` | Long-form STR x patient table |
| `per_str_case_control.tsv` | Descriptive stats per STR: case/control counts, allele means/medians |

**Run**: `qsub compare_gwas_rna.pbs`

---

## 6.3.2.2 — Per-Intervention DEG x STR Crossing

### Step 1: `1_cross_intervention_STRs.py`

Crosses DEGs from each intervention (TSV per comparison) with the cohort STR catalog. Intervention name extracted from filename.

**Inputs**:
- `--deg-dir` — root with `GSE*/` subdirectories (each GSE has intervention TSVs)
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--out-dir` — output directory

**Outputs** (in `--out-dir`):
| File | Description |
|---|---|
| `intervention_strs.tsv` | All STRs annotated by intervention |
| `intervention_outliers.tsv` | STRs with DBSCAN global outliers only |
| `intervention_summary.tsv` | Summary per (intervention, gene): STRs, outliers |

**Run**: `qsub 1_cross_intervention_STRs.pbs`

### Step 2: `2_plot_intervention_summary.R`

Raincloud plot + publication-ready table per intervention.

**Inputs**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--intervention-strs` — `intervention_strs.tsv`
- `--intervention-outliers` — `intervention_outliers.tsv`
- `--intervention-summary` — `intervention_summary.tsv`
- `--out-dir` — output directory

**Outputs** (in `--out-dir`):
| File | Description |
|---|---|
| `intervention_raincloud.png` | Raincloud: allele2_est by intervention, colored case/control |
| `intervention_publication_table.tsv` | Table per intervention |
| `intervention_publication_table.html` | `gt` publication-ready table |

**Run**: `qsub 2_plot_intervention_summary.pbs`

### Step 3: `3_plot_raincloud_per_locus.R`

Raincloud per locus (STR) for outlier DBSCAN samples only. Generates per-GSE and aggregated views.

**Inputs**:
- `--intervention-outliers` — `intervention_outliers.tsv`
- `--out-dir` — output directory

**Outputs** (in `--out-dir`):
| File | Description |
|---|---|
| `ALL/<intervention>_raincloud_per_locus_pXX.png` | Facet across all GSEs |
| `ALL/<intervention>_patients.csv` | Patient-level data |
| `<GSE>/<intervention>_*.png` | Per-GSE views |

**Run**: `qsub 3_submit_raincloud.pbs`

---

## Execution Order

```bash
git pull
qsub 1_cross_intervention_STRs.pbs
# wait for completion
qsub 2_plot_intervention_summary.pbs
# wait for completion
qsub 3_submit_raincloud.pbs
```
