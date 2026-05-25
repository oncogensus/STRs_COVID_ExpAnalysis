cat("\n========================================\n")
cat("LOADING LIBRARIES\n")
cat("========================================\n")

# data.table: High-performance manipulation of large datasets
# rstatix: Pipe-friendly framework for basic statistical tests
# dplyr: Data manipulation grammar (filtering, selecting, mutating)
library(data.table)
library(rstatix)
library(dplyr)

# =========================
# 0. INPUT / OUTPUT PATHS
# =========================

input_file <- "../../samples/STRs_analysis_dataset.tsv"
output_dir <- "results/categorical_data"

# Check if the output directory exists; if not, create it recursively
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
}

# =========================
# 1. LOAD DATA
# =========================

# Fast reading of the TSV file using data.table's fread
merged_dt <- fread(input_file)

# ---- RENAME COLUMNS TO MATCH PIPELINE ----
# Standardize column names for consistency throughout the analysis
setnames(
  merged_dt,
  old = c("sample_id_clean", "pop"),
  new = c("sample_id", "pop_ancestry"),
)

# Calculate the raw total of STRs per sample before any filtering
total_strs_per_sample <- merged_dt[, .(
  total_STRs_raw = .N
), by = sample_id]

# Log dimensions of the loaded dataset
cat("Rows:", nrow(merged_dt), "\n")
cat("Columns:", ncol(merged_dt), "\n")

# =========================
# 2. DERIVE VARIABLES
# =========================

# Calculate the arithmetic mean between the two estimated alleles
merged_dt[, allele_mean := (allele1_est + allele2_est) / 2]

# =========================
# 3. DBSCAN QUALITY FILTER
# =========================

# Define pass criteria: cluster presence, at least one outlier, and low noise ratio
merged_dt[, dbscan_pass := (
  n_clusters > 0 &
  n_outliers >= 1 &
  noise_ratio <= 0.10
)]

cat("STRs passing DBSCAN QC:", sum(merged_dt$dbscan_pass, na.rm = TRUE), "\n")

# Save Quality Control (QC) flags for reproducibility and traceability
fwrite(
  merged_dt[, .(STRs_ID, sample_id, pop_ancestry, dbscan_pass)],
  file.path(output_dir, "dbscan_qc_flags.csv")
)

# Apply filters: keep only records passing QC and with valid ancestry information
filtered_dt <- merged_dt[dbscan_pass == TRUE]
analysis_dt <- filtered_dt[!is.na(pop_ancestry)]

# =========================
# 4. SUMMARY (ALLELES)
# =========================

# Generate descriptive statistics (mean, median, SD) grouped by population ancestry
dist_summary <- analysis_dt[, .(
  mean_a1 = mean(allele1_est, na.rm = TRUE),
  median_a1 = median(allele1_est, na.rm = TRUE),
  sd_a1 = sd(allele1_est, na.rm = TRUE),

  mean_a2 = mean(allele2_est, na.rm = TRUE),
  median_a2 = median(allele2_est, na.rm = TRUE),
  sd_a2 = sd(allele2_est, na.rm = TRUE),

  mean_mean = mean(allele_mean, na.rm = TRUE),
  median_mean = median(allele_mean, na.rm = TRUE),
  sd_mean = sd(allele_mean, na.rm = TRUE)
), by = pop_ancestry]

# =========================
# 5. KRUSKAL-WALLIS (ALLELES)
# =========================

# Perform non-parametric Kruskal-Wallis test to compare distributions across groups
kw_a1 <- kruskal.test(allele1_est ~ pop_ancestry, data = analysis_dt)
kw_a2 <- kruskal.test(allele2_est ~ pop_ancestry, data = analysis_dt)
kw_mean <- kruskal.test(allele_mean ~ pop_ancestry, data = analysis_dt)

# Consolidate statistics and p-values into a results table
kw_results <- data.table(
  variable = c("allele1", "allele2", "allele_mean"),
  statistic = c(kw_a1$statistic, kw_a2$statistic, kw_mean$statistic),
  p_value = c(kw_a1$p.value, kw_a2$p.value, kw_mean$p.value)
)

# =========================
# 6. DUNN POST-HOC (ALLELES)
# =========================

# If Kruskal-Wallis is significant (p < 0.05), run Dunn test to identify group differences
dunn_alleles_list <- list()

if (kw_a1$p.value < 0.05) {
  d <- dunn_test(analysis_dt, allele1_est ~ pop_ancestry, p.adjust.method = "BH")
  d$variable <- "allele1"
  dunn_alleles_list[["a1"]] <- d
}

if (kw_a2$p.value < 0.05) {
  d <- dunn_test(analysis_dt, allele2_est ~ pop_ancestry, p.adjust.method = "BH")
  d$variable <- "allele2"
  dunn_alleles_list[["a2"]] <- d
}

if (kw_mean$p.value < 0.05) {
  d <- dunn_test(analysis_dt, allele_mean ~ pop_ancestry, p.adjust.method = "BH")
  d$variable <- "allele_mean"
  dunn_alleles_list[["mean"]] <- d
}

# Bind all post-hoc results into a single data.table
dunn_alleles_results <- rbindlist(dunn_alleles_list, fill = TRUE)

# =========================
# 7. FORMAT & SAVE ALLELE RESULTS
# =========================

# Pivot data to "long" format for easier visualization (e.g., with ggplot2)
allele_long <- melt(
  analysis_dt,
  id.vars = c("sample_id", "pop_ancestry"),
  measure.vars = c("allele1_est", "allele2_est", "allele_mean"),
  variable.name = "allele_type",
  value.name = "value"
)

# Export statistical summaries and plot-ready data to CSV
fwrite(dist_summary, file.path(output_dir, "alleles_distribution_summary.csv"))
fwrite(kw_results, file.path(output_dir, "alleles_kruskal_results.csv"))
fwrite(allele_long, file.path(output_dir, "plotdata_alleles_long.csv"))
fwrite(analysis_dt, file.path(output_dir, "plotdata_alleles_wide.csv"))

if (nrow(dunn_alleles_results) > 0) {
  fwrite(dunn_alleles_results, file.path(output_dir, "alleles_dunn_results.csv"))
}

# =========================
# 8. DBSCAN OUTLIER BURDEN
# =========================

# Specific QC filter for DBSCAN metrics
merged_dt[, qc_pass := (n_clusters > 0 & noise_ratio <= 0.10)]

# Sample-level aggregation: mantendo os nomes originais para evitar erro nos testes estatísticos
sample_outlier_dt <- merged_dt[qc_pass == TRUE, .(
  total_STRs_after_QC = .N,
  n_outlier_STRs = sum(n_outliers > 0, na.rm = TRUE),
  prop_outliers = mean(n_outliers > 0, na.rm = TRUE), # Removido o sufixo _after_QC
  mean_outlier_strength = mean(abs(as.numeric(outlier_residuals[n_outliers > 0])), na.rm = TRUE) # Removido o sufixo
), by = .(sample_id, pop_ancestry)]

# Merge with the initial raw total counts per sample
sample_outlier_dt <- merge(
  sample_outlier_dt,
  total_strs_per_sample,
  by = "sample_id",
  all.x = TRUE
)

# Data Cleaning: Agora as colunas prop_outliers e mean_outlier_strength existem!
sample_outlier_dt[is.na(mean_outlier_strength), mean_outlier_strength := 0]
sample_outlier_dt <- sample_outlier_dt[!is.na(pop_ancestry)]

# Statistical summary
dist_outlier <- sample_outlier_dt[, .(
  mean_prop = mean(prop_outliers, na.rm = TRUE),
  mean_strength = mean(mean_outlier_strength, na.rm = TRUE),
  n_samples = .N
), by = pop_ancestry]

# =========================
# 9. KRUSKAL-WALLIS (DBSCAN)
# =========================

# Test if outlier proportion and strength significantly differ by ancestry
kw_prop <- kruskal.test(prop_outliers ~ pop_ancestry, data = sample_outlier_dt)
kw_strength <- kruskal.test(mean_outlier_strength ~ pop_ancestry, data = sample_outlier_dt)

kw_dbscan <- data.table(
  variable = c("prop_outliers", "mean_outlier_strength"),
  statistic = c(kw_prop$statistic, kw_strength$statistic),
  p_value = c(kw_prop$p.value, kw_strength$p.value)
)

# =========================
# 10. DUNN POST-HOC (DBSCAN)
# =========================

dunn_dbscan_list <- list()

if (kw_prop$p.value < 0.05) {
  cat("Significant differences in outlier proportion. Running Dunn test...\n")
  d_prop <- dunn_test(sample_outlier_dt, prop_outliers ~ pop_ancestry, p.adjust.method = "BH")
  d_prop$variable <- "prop_outliers"
  dunn_dbscan_list[["prop"]] <- d_prop
}

if (kw_strength$p.value < 0.05) {
  cat("Significant differences in outlier strength. Running Dunn test...\n")
  d_strength <- dunn_test(sample_outlier_dt, mean_outlier_strength ~ pop_ancestry, p.adjust.method = "BH")
  d_strength$variable <- "mean_outlier_strength"
  dunn_dbscan_list[["strength"]] <- d_strength
}

dunn_dbscan_results <- rbindlist(dunn_dbscan_list, fill = TRUE)

# =========================
# 11. SAVE DBSCAN RESULTS
# =========================

# Save DBSCAN-related summaries and statistical results
fwrite(dist_outlier, file.path(output_dir, "dbscan_distribution_summary.csv"))
fwrite(kw_dbscan, file.path(output_dir, "dbscan_kruskal_results.csv"))

if (nrow(dunn_dbscan_results) > 0) {
  fwrite(dunn_dbscan_results, file.path(output_dir, "dbscan_dunn_results.csv"))
}

# Prepare DBSCAN metric data in long format for plotting
dbscan_long <- melt(
  sample_outlier_dt,
  id.vars = c("sample_id", "pop_ancestry"),
  measure.vars = c("prop_outliers", "mean_outlier_strength"),
  variable.name = "metric",
  value.name = "value"
)

fwrite(dbscan_long, file.path(output_dir, "plotdata_dbscan_long.csv"))
fwrite(sample_outlier_dt, file.path(output_dir, "plotdata_dbscan_wide.csv"))

cat("\n========================================\n")
cat("ANALYSIS COMPLETED: DATA SAVED FOR PLOTTING\n")
cat("========================================\n")