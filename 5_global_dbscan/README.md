# 5 — Global DBSCAN Analysis

Normalizes the largest STR allele by covariates and detects outlier samples per locus with DBSCAN.

## Purpose

Remove the confounding effect of sex, ancestry (PCA), sequencing depth and age from `allele2_est`, then flag outlier STR alleles genome-wide using density-based clustering.

## Structure

```
5_global_dbscan/
├── norm_test/
│   └── 5.1_norm_dbscan.r      # Regression-based normalization
└── outliers_search/
    └── 5.2_dbscan_str.r       # Per-locus DBSCAN outlier detection
```

## 5.1 — STR Normalization

Regresses `allele2_est ~ sex + EV1 + EV2 + EV3 + depth + age`, stratified per `STRs_ID`, and stores the residuals.

- Regression is only computed for loci with **n > 10** valid samples; otherwise residuals are set to `NA`.
- IDs are standardized (`normalize_ids`): BAM suffixes removed, replicate suffixes dropped, leading zeros normalized.

### Execution

```bash
cd 5_global_dbscan/norm_test
Rscript 5.1_norm_dbscan.r
```

### Outputs

- `norm_test/STRs_normalized_residuals.tsv` — residuals per STR/sample

## 5.2 — DBSCAN Outlier Detection

Applies DBSCAN on the residuals of each STR independently:

- `minPts = max(2, ceil(log2(2 * n_valid)))`
- `eps` derived from the residual spread (5th/95th percentile range), floored at `1e-6`
- Cutoff = largest residual among non-noise points (minimum 2); loci with a single cluster or no noise get `cutoff = Inf`

### Execution

```bash
cd 5_global_dbscan/outliers_search
Rscript 5.2_dbscan_str.r
```

### Outputs

- `results_dbscan/outliers_per_str.tsv` — one row per STR: sample/residual lists, `n_clusters`, `noise_ratio`

## Environment

- `dbscan-r` (micromamba): `data.table`, `dbscan` (pinned as `dbscan-r.yaml`, repo root)

**Reference**: DBSCAN approach based on [AD_STR](https://github.com/mhguo1/AD_STR/tree/main).