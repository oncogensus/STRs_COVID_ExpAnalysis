#!/usr/bin/env Rscript
# 1_burden_analysis.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Complementary burden analyses for reviewer comment 1.10 (association
#   model). Compares the RELATIVE burden of DBSCAN outlier STRs per individual
#   between fatal COVID-19 cases and survivors:
#     (1) Global relative burden                -> Mann-Whitney U
#     (2) Global relative burden (adjusted)     -> Firth logistic regression
#                                                  fatal ~ burden + age + sex + PC1
#     (3) Relative burden within DEGs           -> Mann-Whitney U
#     (4) Relative burden per genomic region    -> Mann-Whitney U
#   All analyses are exploratory; p-values are nominal (no FDR).
#
# INPUTS (command-line arguments; defaults are repo-relative)
#   --str-catalog   samples/STRs_analysis_dataset.tsv (stage 6.1)
#   --deg-strs      intervention_strs.tsv (stage 6.3.2.2 step 1)
#   --pca           EthSEQ_Results_3D/Report.PCAcoord (stage 4)
#   --out-dir       Output directory (default: results)
#
# OUTPUTS (in --out-dir)
#   burden_per_sample.csv     Per-individual absolute/relative burden
#   burden_global_mw.csv      Panel A: global relative burden (Mann-Whitney)
#   burden_firth.csv          Panel B: Firth logistic regression
#   burden_deg_mw.csv         Panel C: relative burden within DEGs
#   burden_region_mw.csv      Panel D: relative burden per genomic region
#
# ENVIRONMENT
#   r_enrich_env (micromamba): data.table, dplyr, stringr, logistf
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# logistf is only required for the Firth logistic regression (Panel B).
has_logistf <- requireNamespace("logistf", quietly = TRUE)
if (!has_logistf) {
  cat("[WARN] Package 'logistf' not available. Panel B (Firth regression) will be skipped.\n")
  cat("[WARN] Install with: install.packages('logistf') or micromamba install -n r_enrich_env -c conda-forge r-logistf\n")
}

# ==========================================
# Parse arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  return(default)
}

path_str_catalog <- parse_arg("--str-catalog", "../../../samples/STRs_analysis_dataset.tsv")
path_deg_strs    <- parse_arg("--deg-strs",    "../6.3.2_RNA_data_analysis/6.3.2.2_RNA_matrix/results/intervention_strs.tsv")
path_pca         <- parse_arg("--pca",         "../../../4_ancestry/EthSEQ_Results_3D/Report.PCAcoord")
out_dir          <- parse_arg("--out-dir",     "results")

if (is.null(path_str_catalog)) stop("Missing argument: --str-catalog")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# Constants
# ==========================================
# Outcome labels mapped to fatal (1) vs survivor (0). Cohort uses case/control.
FATAL_LABELS    <- c("case", "fatal", "death", "obito", "óbito", "nonsurvivor", "non-survivor")
SURVIVOR_LABELS <- c("control", "survivor", "surv", "sobrevivente")

# Stored region value -> manuscript label.
REGION_LABELS <- c(
  intergenic       = "Intergenic",
  intron           = "Intron",
  others           = "Non-Coding Elements",
  promoter         = "Promoter",
  three_prime_utr  = "3' UTR",
  non_coding_exons = "Non-coding Exons",
  five_prime_utr   = "5' UTR",
  CDS              = "CDS"
)

# ==========================================
# Helpers
# ==========================================
normalize_ids <- function(ids) {
  ids %>%
    as.character() %>%
    toupper() %>%
    str_remove("(?i)[._-]?\\d*BAM.*$") %>%
    str_remove("(?<=\\d)-[0-9]$") %>%
    str_replace_all("-0+([0-9]+)", "-\\1") %>%
    str_trim()
}

# Format a numeric vector as "median [Q1-Q3]".
med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.4f [%.4f-%.4f]", q[2], q[1], q[3])
}

# Mann-Whitney U of metric_col between fatal (1) and survivor (0).
mw_summary <- function(dt, metric_col, group_col = "outcome") {
  x <- dt[get(group_col) == 1][[metric_col]]
  y <- dt[get(group_col) == 0][[metric_col]]
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 1 || length(y) < 1) {
    return(data.table(n_survivors = length(y), n_fatal = length(x),
                      survivors_median_IQR = NA_character_,
                      fatal_median_IQR = NA_character_, p = NA_real_))
  }
  p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  data.table(
    n_survivors          = length(y),
    n_fatal              = length(x),
    survivors_median_IQR = med_iqr(y),
    fatal_median_IQR     = med_iqr(x),
    p                    = p
  )
}

# ==========================================
# 1. Load data
# ==========================================
cat("--- Burden analyses: STR outliers x COVID-19 fatality ---\n")
cat(sprintf("  str-catalog: %s\n", path_str_catalog))
merged_dt <- fread(path_str_catalog)

# Strip the DBSCAN global suffix (as in 6.5) for shorter column names.
db_cols <- grep("_dbscan_global$", names(merged_dt), value = TRUE)
if (length(db_cols) > 0) {
  setnames(merged_dt, db_cols, sub("_dbscan_global$", "", db_cols))
}
cat(sprintf("  Rows: %d | Columns: %d\n", nrow(merged_dt), ncol(merged_dt)))

# ==========================================
# 2. DBSCAN QC + outcome coding
# ==========================================
# Same QC criterion used in the ancestry outlier-burden section (6.5).
merged_dt[, qc_pass := (
  !is.na(n_clusters) & n_clusters > 0 &
  !is.na(noise_ratio) & noise_ratio <= 0.10
)]
merged_dt[, is_outlier := qc_pass & !is.na(n_outliers) & n_outliers > 0]
merged_dt[is.na(is_outlier), is_outlier := FALSE]

g <- tolower(trimws(as.character(merged_dt$group)))
g_num <- suppressWarnings(as.integer(g))
merged_dt[, outcome := NA_integer_]
merged_dt[!is.na(g_num) & g_num %in% c(0L, 1L), outcome := g_num]
merged_dt[is.na(outcome) & g %in% FATAL_LABELS, outcome := 1L]
merged_dt[is.na(outcome) & g %in% SURVIVOR_LABELS, outcome := 0L]

cat(sprintf("  Outcome labels observed: %s\n",
            paste(sort(unique(as.character(merged_dt$group))), collapse = ", ")))
cat(sprintf("  Samples with outcome: %d\n",
            length(unique(merged_dt[!is.na(outcome)]$sample_id))))
if (any(is.na(merged_dt$outcome))) {
  cat(sprintf("  [WARN] Unmapped outcome labels: %s\n",
              paste(sort(unique(as.character(merged_dt[is.na(outcome)]$group))), collapse = ", ")))
}

# ==========================================
# 3. Per-sample global relative burden
# ==========================================
sample_dt <- merged_dt[, .(
  total_STRs     = sum(qc_pass, na.rm = TRUE),
  n_outlier_STRs = sum(is_outlier, na.rm = TRUE)
), by = sample_id]
sample_dt[, burden_rel := n_outlier_STRs / pmax(total_STRs, 1)]

# One phenotype record per sample.
pheno <- merged_dt[, .(outcome = outcome[1], age = age[1], sex = sex[1]), by = sample_id]
pheno <- pheno[!is.na(outcome)]

# PC1 from EthSEQ ancestry coordinates.
pca <- fread(path_pca)
setnames(pca, names(pca)[1], "sample_id_clean")
pca[, sample_id_clean := normalize_ids(sample_id_clean)]
pca <- pca[, .(sample_id_clean, PC1 = EV1)]

burden <- merge(sample_dt, pheno, by = "sample_id")
burden[, sample_id_clean := normalize_ids(sample_id)]
burden <- merge(burden, pca, by = "sample_id_clean", all.x = TRUE)
cat(sprintf("  Samples in burden table: %d | with PC1: %d\n",
            nrow(burden), sum(!is.na(burden$PC1))))

# ==========================================
# 4. Analysis 1: global burden (Mann-Whitney)
# ==========================================
cat("\n[1] Global relative burden - Mann-Whitney U\n")
global_mw <- mw_summary(burden, "burden_rel")
global_mw[, context := "Global"]
print(global_mw)

# ==========================================
# 5. Analysis 2: Firth logistic regression
# ==========================================
cat("\n[2] Firth logistic regression: fatal ~ burden_rel + age + sex + PC1\n")
fit_res <- NULL
if (has_logistf) {
  fit_dt <- as.data.frame(burden[!is.na(outcome) & !is.na(age) & !is.na(sex) & !is.na(PC1)])
  fit_dt$sex <- as.factor(fit_dt$sex)
  if (nrow(fit_dt) > 0 && length(unique(fit_dt$outcome)) == 2) {
    fit <- tryCatch(
      logistf::logistf(outcome ~ burden_rel + age + sex + PC1, data = fit_dt),
      error = function(e) { cat("[WARN] Firth model did not converge:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(fit)) {
      fit_res <- data.table(
        predictor = names(coef(fit)),
        OR        = exp(coef(fit)),
        CI_low    = exp(fit$ci.lower),
        CI_high   = exp(fit$ci.upper),
        p         = fit$prob
      )
      print(fit_res)
    } else {
      cat("Not estimable: model did not converge.\n")
    }
  } else {
    cat("[WARN] Insufficient data for Firth model (need both outcomes and complete covariates).\n")
  }
}

# ==========================================
# 6. Analysis 3: relative burden within DEGs
# ==========================================
cat("\n[3] Relative burden within DEGs - Mann-Whitney U\n")
deg_mw <- NULL
if (file.exists(path_deg_strs)) {
  deg_strs <- unique(fread(path_deg_strs)$STRs_ID)
  cat(sprintf("  STRs within DEGs: %d\n", length(deg_strs)))
  merged_dt[, is_deg := STRs_ID %in% deg_strs]
  deg_dt <- merged_dt[qc_pass == TRUE & is_deg == TRUE, .(
    total_deg = .N,
    n_out_deg = sum(is_outlier, na.rm = TRUE)
  ), by = sample_id]
  deg_dt[, burden_rel_DEG := n_out_deg / pmax(total_deg, 1)]
  deg_burden <- merge(deg_dt, pheno, by = "sample_id")
  deg_mw <- mw_summary(deg_burden, "burden_rel_DEG")
  deg_mw[, context := "DEGs"]
  print(deg_mw)
} else {
  cat(sprintf("  [WARN] DEG file not found: %s - Panel C skipped.\n", path_deg_strs))
}

# ==========================================
# 7. Analysis 4: relative burden by genomic region
# ==========================================
cat("\n[4] Relative burden by genomic region - Mann-Whitney U\n")
regions <- sort(unique(merged_dt[qc_pass == TRUE]$region))
region_list <- lapply(regions, function(r) {
  sub <- merged_dt[qc_pass == TRUE & region == r]
  agg <- sub[, .(
    total_reg = .N,
    n_out_reg = sum(is_outlier, na.rm = TRUE)
  ), by = sample_id]
  agg[, burden_rel_reg := n_out_reg / pmax(total_reg, 1)]
  agg <- merge(agg, pheno, by = "sample_id")
  mw <- mw_summary(agg, "burden_rel_reg")
  lbl <- unname(REGION_LABELS[r])
  if (is.na(lbl)) lbl <- r
  mw[, region := lbl]
  mw[, not_estimable := (sum(sub$is_outlier, na.rm = TRUE) == 0)]
  mw
})
region_mw <- rbindlist(region_list)
setcolorder(region_mw, c("region", "n_survivors", "n_fatal",
                         "survivors_median_IQR", "fatal_median_IQR", "p", "not_estimable"))
print(region_mw)

# ==========================================
# 8. Save outputs
# ==========================================
cat("\nSaving results...\n")
fwrite(burden,    file.path(out_dir, "burden_per_sample.csv"))
fwrite(global_mw, file.path(out_dir, "burden_global_mw.csv"))
if (!is.null(fit_res)) {
  fwrite(fit_res, file.path(out_dir, "burden_firth.csv"))
}
if (!is.null(deg_mw)) {
  fwrite(deg_mw, file.path(out_dir, "burden_deg_mw.csv"))
}
fwrite(region_mw, file.path(out_dir, "burden_region_mw.csv"))

cat("\nDone.\n")
