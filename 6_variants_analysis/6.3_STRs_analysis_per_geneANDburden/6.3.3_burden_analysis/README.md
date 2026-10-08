# 6.3.3 — Burden Analysis

Complementary burden analyses for the association-model reviewer comment (1.10). Compares the **relative burden** of DBSCAN outlier STRs per individual between fatal COVID-19 cases and survivors.

## Objective

Provide a parsimonious, exploratory set of burden tests that complement the locus-specific Mann-Whitney U tests (6.3.2 step 4). P-values are reported as nominal. A Benjamini–Hochberg FDR (`p_adj`) is computed for panels B (across the non-intercept predictors), C (across DEG contexts) and D (across estimable regions); panel A is a single comparison and is shown with nominal *P* only. Four analyses are performed:

| # | Context | Test |
|---|---|---|
| 1 | Global relative burden | Mann-Whitney U |
| 2 | Global relative burden (adjusted) | Firth logistic regression: `fatal ~ burden_pct + age + sex + PC1` (burden per 1%, PC1 per SD) |
| 3 | Relative burden within DEGs (global + per intervention) | Mann-Whitney U |
| 4 | Relative burden per genomic region | Mann-Whitney U |

**Definitions**

- Relative burden = `n_outliers / n_total_strs` per individual (QC: `n_clusters > 0`, `noise_ratio <= 0.10`).
- Relative burden within DEGs = outliers restricted to STRs within DEGs / total STRs within DEGs. Reported globally (all DEGs) and per intervention.
- Relative burden per region = outliers in the region / total STRs in the region (8 regions: Intergenic, Intron, Non-Coding Elements, Promoter, 3′ UTR, Non-coding Exons, 5′ UTR, CDS).
- Outcome: fatal (case) = 1 vs survivor (control) = 0.

## Structure

```
6.3.3_burden_analysis/
├── 1_burden_analysis.R      # Burden analyses (4 tests)
└── 2_burden_tables.R        # Supplementary tables (gt HTML)
```

## Inputs

| File | Source | Description |
|---|---|---|
| `STRs_analysis_dataset.tsv` | `samples/` (stage 6.1) | Unified STR dataset (group, age, sex, DBSCAN global metrics) |
| `intervention_strs.tsv` | `6.3.2.2_RNA_matrix/results/` | STRs within DEGs (per intervention), for Panel C |
| `Report.PCAcoord` | `4_ancestry/EthSEQ_Results_3D/` | Ancestry coordinates (EV1 = PC1) |

## Outputs (`results/`)

| File | Description |
|---|---|
| `burden_per_sample.csv` | Per-individual absolute and relative burden |
| `burden_global_mw.csv` | Panel A: global relative burden (Mann-Whitney); includes BH FDR (`p_adj`) |
| `burden_firth.csv` | Panel B: Firth logistic regression (OR, 95% CI, p, `p_adj`) |
| `burden_deg_mw.csv` | Panel C: relative burden within DEGs (global + per intervention); includes BH FDR (`p_adj`) |
| `burden_region_mw.csv` | Panel D: relative burden per genomic region; includes BH FDR (`p_adj`) |
| `tables/burden_analysis.html` | All four panels (gt HTML) |
| `tables/burden_global.html`, `burden_firth.html`, `burden_deg.html`, `burden_region.html` | Individual panels (gt HTML) |
| `tables/legends.md` | Table legends (external; not embedded in the HTML) |

Tables omit comparisons without variance in both groups (IQR span = 0); panel A shows nominal *P* only (no FDR). Legends are written to `legends.md` for manual placement in the manuscript.

## Execution

```bash
cd 6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.3.3_burden_analysis
qsub 1_burden_analysis.pbs          # on the cluster (r_enrich_env)
qsub 2_burden_tables.pbs            # supplementary tables (gt HTML)
# or directly:
Rscript 1_burden_analysis.R
Rscript 2_burden_tables.R --results-dir results --out-dir results/tables
```

Paths are repo-relative by default and can be overridden with `--str-catalog`, `--deg-strs`, `--pca`, and `--out-dir`. Requires `intervention_strs.tsv` (stage 6.3.2.2 step 1) for Panel C; if missing, Panel C is skipped with a warning.

## Environment

- `r_enrich_env` (micromamba): `data.table`, `dplyr`, `stringr`, `logistf` (pinned as `r_enrich_env.yaml`, repo root)

> Panel B (Firth) requires `logistf`. If unavailable, the script skips Panel B with a warning and the remaining analyses still run.
