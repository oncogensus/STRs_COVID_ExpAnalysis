#!/usr/bin/env Rscript
# 7_outlier_report.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Relatorio unificado de outliers DBSCAN para a coorte global:
#     - General Summary (total, signal, no_signal)
#     - Cluster Distribution (1, 2, 3+ clusters)
#     - Noise Tiers (0-5%, 5-10%, 10-20%, 20-50%, >50%)
#     - Outlier Frequency (0 vs >0 outliers)
#   Inclui tabela publication-ready (gt HTML).
#   Avalia TODOS os registros (sem filtro de outliers).
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog    STRs_analysis_dataset.tsv (coorte global)
#   --out-dir        Diretorio de saida
#
# SAIDAS
#   unified_binary_outlier_report.csv
#   unified_technical_report.csv
#   Technical_Validation_Table.html
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(gt)
})

# ==========================================
# Parse arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  return(default)
}

path_str_catalog <- parse_arg("--str-catalog")
out_dir          <- parse_arg("--out-dir", ".")

if (is.null(path_str_catalog)) {
  stop("Argumento ausente: --str-catalog")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("--- Unified Outlier Report (Cohort Global) ---\n")

df_strs <- fread(path_str_catalog, header = TRUE, sep = "\t")
cat(sprintf("  STRs_analysis_dataset.tsv: %d linhas\n", nrow(df_strs)))

# Rename DBSCAN columns (with _dbscan_global suffix) to short names
col_names <- names(df_strs)
if ("n_clusters_dbscan_global" %in% col_names) {
  setnames(df_strs,
           c("n_clusters_dbscan_global", "noise_ratio_dbscan_global",
             "n_outliers_dbscan_global", "outlier_samples_dbscan_global",
             "outlier_residuals_dbscan_global"),
           c("n_clusters", "noise_ratio",
             "n_outliers", "outlier_samples",
             "outlier_residuals"))
  cat("  Colunas DBSCAN renomeadas (removido _dbscan_global)\n")
}

# ==========================================
# 2. Build unified report
# ==========================================
cat("\nGerando relatorio unificado...\n")

# A. General Summary
general_summary <- df_strs[, .(
  total_obs   = .N,
  with_signal = sum(!is.na(n_clusters) & n_clusters > 0),
  no_signal   = .N - sum(!is.na(n_clusters) & n_clusters > 0)
)]

df_general <- melt(general_summary, variable.name = "metric", value.name = "raw_value")
df_general[, category := "General Summary"]
df_general[, percentage := round(raw_value / first(raw_value) * 100, 2)]

# B. Cluster Distribution (signal only)
signal_data <- df_strs[!is.na(n_clusters) & n_clusters > 0]

cluster_dist <- signal_data[, .N, by = n_clusters]
setnames(cluster_dist, "N", "raw_value")
cluster_dist[, metric := paste0(n_clusters, " Cluster(s)")]
cluster_dist[, category := "Cluster Distribution"]
cluster_dist[, percentage := round(raw_value / sum(raw_value) * 100, 2)]
cluster_dist <- cluster_dist[, .(category, metric, raw_value, percentage)]

# C. Noise Tiers (signal only)
noise_data <- signal_data[, .(
  metric = case_when(
    noise_ratio <= 0.05 ~ "Noise 00-05%",
    noise_ratio <= 0.10 ~ "Noise 05-10%",
    noise_ratio <= 0.20 ~ "Noise 10-20%",
    noise_ratio <= 0.50 ~ "Noise 20-50%",
    TRUE                ~ "Noise > 50%"
  )
)]

noise_tiers <- noise_data[, .N, by = metric]
setnames(noise_tiers, "N", "raw_value")
noise_tiers[, category := "Noise Tiers"]
noise_tiers[, percentage := round(raw_value / sum(raw_value) * 100, 2)]

# D. Outlier Frequency (0 vs >0)
outlier_freq <- signal_data[, .(
  has_outlier = ifelse(n_outliers > 0,
                       "Observations with >0 Outliers",
                       "Observations with 0 Outliers")
)]

outlier_counts <- outlier_freq[, .N, by = has_outlier]
setnames(outlier_counts, c("has_outlier", "N"), c("metric", "raw_value"))
outlier_counts[, category := "Outlier Frequency"]
outlier_counts[, percentage := round(raw_value / sum(raw_value) * 100, 2)]

# Bind all
df_report <- rbind(df_general, cluster_dist, noise_tiers, outlier_counts)
df_report <- df_report[, .(category, metric, raw_value, percentage)]

cat("\n--- Relatorio Unificado ---\n")
print(as.data.frame(df_report))

# ==========================================
# 3. Export CSV
# ==========================================
out_csv <- file.path(out_dir, "unified_binary_outlier_report.csv")
fwrite(df_report, out_csv, sep = ";", dec = ",")
cat(sprintf("\nCSV salvo em: %s\n", out_csv))

# ==========================================
# 4. Publication-ready table (gt HTML)
# ==========================================
cat("\nGerando tabela publication-ready...\n")

# Friendly labels
df_friendly <- copy(df_report)
df_friendly[, metric := fifelse(
  metric == "total_obs",     "Total Observations",
  fifelse(metric == "with_signal",   "Observations with Signal",
  fifelse(metric == "no_signal",     "Observations without Signal",
  fifelse(metric == "1 Cluster(s)",  "1 Cluster", metric))))]
df_friendly[, percentage_str := paste0(round(percentage, 2), "%")]

pub_gt <- df_friendly %>%
  group_by(category) %>%
  gt() %>%
  tab_header(
    title = md("**Table 1. Technical Validation of STR Genotyping**"),
    subtitle = "DBSCAN efficiency, cluster stability, and noise metrics"
  ) %>%
  cols_label(
    metric        = "Metric",
    raw_value     = "Absolute Frequency (n)",
    percentage_str = "Relative Frequency (%)"
  ) %>%
  fmt_number(
    columns = c(raw_value),
    decimals = 0,
    use_seps = TRUE
  ) %>%
  tab_options(
    table.width = px(750),
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.background.color = "#F9F9F9",
    table.font.size = px(14),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.align = "left",
    data_row.padding = px(6)
  ) %>%
  cols_align(align = "left",   columns = c(metric)) %>%
  cols_align(align = "center", columns = c(raw_value, percentage_str)) %>%
  tab_source_note(
    source_note = md("*Note: Signal detection and noise metrics are calculated based on DBSCAN clustering parameters.*")
  )

out_gt_html <- file.path(out_dir, "Technical_Validation_Table.html")
gtsave(pub_gt, out_gt_html)
cat(sprintf("Tabela gt HTML salva em: %s\n", out_gt_html))

# Save raw CSV for the report
out_report_csv <- file.path(out_dir, "unified_technical_report.csv")
fwrite(df_report, out_report_csv, sep = ";", dec = ",")
cat(sprintf("CSV exportado em: %s\n", out_report_csv))

cat("\nConcluido.\n")
