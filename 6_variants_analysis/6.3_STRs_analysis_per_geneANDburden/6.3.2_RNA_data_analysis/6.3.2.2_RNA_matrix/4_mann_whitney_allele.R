#!/usr/bin/env Rscript
# 4_mann_whitney_allele.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Mann-Whitney U test: case vs. control para STRs com outliers DBSCAN,
#   por intervencao (GSE). Dois testes:
#     (a) allele2_est com N >= 3 por grupo
#     (b) mean_allele = (allele1 + allele2) / 2 com N >= 3 por grupo
#   Inclui effect size (rank-biserial r) e FDR (Benjamini-Hochberg).
#   Gera tabela publication-ready (gt HTML) com todos os loci testados.
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog    intervention_outliers.tsv (saida do step 1)
#   --intervention   GSE (ex: GSE157103) ou ALL (default: ALL)
#   --out-dir        Diretorio de saida
#
# SAIDAS
#   {intervention}_allele2_mw.csv
#   {intervention}_mean_allele_mw.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rstatix)
  library(gt)
  library(scales)
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
intervention     <- parse_arg("--intervention", "ALL")
out_dir          <- parse_arg("--out-dir", ".")

if (is.null(path_str_catalog)) {
  stop("Argumento ausente: --str-catalog")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("--- Mann-Whitney U Test: Case vs Control (Outlier Loci) ---\n")

df <- fread(path_str_catalog, header = TRUE, sep = "\t")
setnames(df, "intervention", "comparison_type")
cat(sprintf("  intervention_outliers.tsv: %d linhas\n", nrow(df)))

# ==========================================
# 2. Filter by DBSCAN quality (stricter than file default)
# ==========================================
df_filtered <- df[
  n_clusters_dbscan_global >= 1 &
  noise_ratio_dbscan_global < 0.10 &
  n_outliers_dbscan_global >= 1
]
cat(sprintf("  Apos filtro DBSCAN (clusters >= 1, noise < 10%%, outliers >= 1): %d linhas\n",
            nrow(df_filtered)))

# ==========================================
# 3. Filter by intervention
# ==========================================
if (intervention != "ALL") {
  df_filtered <- df_filtered[gse == intervention]
  cat(sprintf("  Intervention '%s': %d linhas\n", intervention, nrow(df_filtered)))
} else {
  cat(sprintf("  ALL interventions: %d linhas\n", nrow(df_filtered)))
}

if (nrow(df_filtered) == 0) {
  cat("Nenhuma observacao apos filtragem. Saindo.\n")
  quit(status = 0)
}

# ==========================================
# 4. Ensure numeric types
# ==========================================
df_filtered[, allele2_est := as.numeric(allele2_est)]
df_filtered[, allele1_est := as.numeric(allele1_est)]
df_filtered[, group := factor(group, levels = c("case", "control"))]

# ==========================================
# Helper: run Mann-Whitney for a given metric
# ==========================================
run_mann_whitney <- function(data, metric_col, min_n_per_group, label) {
  cat(sprintf("\n>>> %s (N per group >= %d) <<<\n", label, min_n_per_group))

  # Prepare metric
  dt <- copy(data)
  dt[, target_metric := get(metric_col)]
  dt <- dt[!is.na(target_metric)]

  cat(sprintf("  Observacoes validas: %d\n", nrow(dt)))

  # Filter loci with sufficient samples per group (per variant x study)
  eligible <- dt[, .N, by = .(STRs_ID, gse, group)]
  eligible <- dcast(eligible, STRs_ID + gse ~ group, value.var = "N", fill = 0)
  eligible <- eligible[case >= min_n_per_group & control >= min_n_per_group]

  cat(sprintf("  Combinacoes variantexestudo elegiveis: %d\n", nrow(eligible)))

  if (nrow(eligible) == 0) {
    cat("  Nenhuma combinacao elegivel.\n")
    return(data.table())
  }

  # Run Wilcoxon test per variant x study
  dt_eligible <- merge(dt, eligible[, .(STRs_ID, gse)],
                       by = c("STRs_ID", "gse"), allow.cartesian = TRUE)

  results <- dt_eligible %>%
    group_by(STRs_ID, gse) %>%
    wilcox_test(target_metric ~ group) %>%
    ungroup() %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    as.data.table()

  # Add sample counts per group
  n_counts <- dt_eligible[, .(n_total = .N), by = .(STRs_ID, gse, group)]
  n_wide <- dcast(n_counts, STRs_ID + gse ~ group, value.var = "n_total", fill = 0)
  setnames(n_wide, c("case", "control"), c("n_cases", "n_controls"))
  results <- merge(results, n_wide, by = c("STRs_ID", "gse"), all.x = TRUE)

  # Compute effect size manually: rank-biserial r = 1 - (2*U) / (n1*n2)
  results[, effsize := 1 - (2 * statistic) / (n_cases * n_controls)]

  results <- results[order(p.adj, STRs_ID, gse)]

  cat(sprintf("  Total de testes: %d\n", nrow(results)))
  cat(sprintf("  Significativos (p.adj < 0.05): %d\n", sum(results$p.adj < 0.05, na.rm = TRUE)))

  return(results)
}

# ==========================================
# 5. Run tests
# ==========================================
# (a) Allele 2
res_allele2 <- run_mann_whitney(df_filtered, "allele2_est", min_n_per_group = 3,
                                label = "Mann-Whitney: Allele 2")

# (b) Mean Allele
df_filtered[, mean_allele := (allele1_est + allele2_est) / 2]
res_mean <- run_mann_whitney(df_filtered, "mean_allele", min_n_per_group = 3,
                             label = "Mann-Whitney: Mean Allele")

# ==========================================
# 6. Export individual CSVs
# ==========================================
suffix <- if (intervention == "ALL") "ALL" else intervention

if (nrow(res_allele2) > 0) {
  out_allele2 <- file.path(out_dir, paste0(suffix, "_allele2_mw.csv"))
  fwrite(as.data.frame(res_allele2), out_allele2, sep = ";", dec = ",")
  cat(sprintf("\n  Allele2 MW salvo em: %s\n", out_allele2))
} else {
  cat("\n  Nenhum resultado para Allele 2.\n")
}

if (nrow(res_mean) > 0) {
  out_mean <- file.path(out_dir, paste0(suffix, "_mean_allele_mw.csv"))
  fwrite(as.data.frame(res_mean), out_mean, sep = ";", dec = ",")
  cat(sprintf("  Mean Allele MW salvo em: %s\n", out_mean))
} else {
  cat("  Nenhum resultado para Mean Allele.\n")
}

# ==========================================
# 7. Combine results and build publication-ready table
# ==========================================
cat("\n--- Publication-Ready Table ---\n")

intervention_labels <- c(
  "GSE183533" = "Fatal COVID-19 vs. Controls",
  "GSE188847" = "Non-Survivors vs. Controls",
  "GSE157103" = "ICU vs. Non-Critical"
)

fmt_pval <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
         ifelse(x < 0.05, sprintf("%.3f", x),
                sprintf("%.2f", x)))
}

# Combine both metrics
res_allele2 <- res_allele2[, Metric := "Allele 2"]
res_mean <- res_mean[, Metric := "Mean Allele"]
res_combined <- rbind(res_allele2, res_mean, fill = TRUE)

# Map gene_name from STRs_ID
gene_map <- unique(df_filtered[, .(STRs_ID, gene_name)])
res_combined <- merge(res_combined, gene_map, by = "STRs_ID", all.x = TRUE)

if (nrow(res_combined) == 0) {
  cat("  Nenhum resultado para tabela.\n")
  cat("\nConcluido.\n")
  quit(status = 0)
}

# Extract motif:size from STRs_ID
res_combined[, Variant := sub("^[^:]+:[^:]+:(.+)$", "\\1", STRs_ID)]

# Map intervention labels
res_combined[, Comparison := fifelse(
  gse %in% names(intervention_labels),
  intervention_labels[gse], gse
)]

# Build table
tbl <- res_combined[, .(
  Comparison = Comparison,
  Gene = gene_name,
  Variant = Variant,
  Metric = Metric,
  `N Cases` = n_cases,
  `N Controls` = n_controls,
  `U Statistic` = round(statistic, 2),
  `p-value` = fmt_pval(p),
  FDR = fmt_pval(p.adj),
  `Effect Size (r)` = round(effsize, 3),
  Sig = fifelse(p.adj < 0.05, "**", fifelse(p < 0.05, "*", ""))
)]
tbl <- tbl[order(Comparison, Gene, Metric)]

n_loci <- uniqueN(tbl$Gene)
n_comps <- uniqueN(tbl$Comparison)

cat(sprintf("  Tabela: %d testes em %d loci x %d comparacoes\n", nrow(tbl), n_loci, n_comps))

# Build gt table
out_gt <- tbl %>%
  gt(groupname_col = "Comparison") %>%
  tab_header(
    title = md("**Mann-Whitney U Test: Case vs Control**"),
    subtitle = sprintf("%d loci tested across %d comparisons", n_loci, n_comps)
  ) %>%
  tab_spanner(
    label = "Sample Size",
    columns = c(`N Cases`, `N Controls`)
  ) %>%
  tab_spanner(
    label = md("Statistical Test^a^"),
    columns = c(`U Statistic`, `p-value`, FDR, Sig)
  ) %>%
  tab_spanner(
    label = md("Effect Size^b^"),
    columns = c(`Effect Size (r)`)
  ) %>%
  cols_label(
    Gene = "Gene",
    Variant = "Variant",
    Metric = "Metric",
    `N Cases` = "Cases",
    `N Controls` = "Controls",
    `U Statistic` = "U",
    `p-value` = "p-value",
    FDR = "FDR",
    `Effect Size (r)` = md("r^b^"),
    Sig = md("Sig.^c^")
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(
      cell_text(weight = "bold", size = px(10)),
      cell_borders(sides = "bottom", weight = px(1.5), color = "grey60")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(10)),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(size = px(9)),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(9)),
    locations = cells_body(columns = Gene)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "grey95"),
      cell_text(weight = "bold", size = px(11))
    ),
    locations = cells_row_groups()
  ) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(9),
    heading.align = "left",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.width = px(1.5),
    table_body.border.bottom.color = "black",
    table_body.hlines.color = "grey92",
    table_body.hlines.width = px(0.5),
    table.border.left.width = px(0),
    table.border.right.width = px(0),
    data_row.padding = px(3),
    row_group.padding = px(6)
  ) %>%
  tab_source_note(
    source_note = md("*^a^* Mann-Whitney U test with Benjamini\u2013Hochberg FDR correction. Groups: case vs control per variant per study.")
  ) %>%
  tab_source_note(
    source_note = md("*^b^* Effect size: rank-biserial correlation (r). Interpretation: |r| < 0.1 negligible, 0.1\u20130.3 small, 0.3\u20130.5 medium, > 0.5 large.")
  ) %>%
  tab_source_note(
    source_note = md("*^c^* Significance: ** p < 0.05 (FDR); * p < 0.05 (nominal, uncorrected).")
  )

out_html <- file.path(out_dir, paste0(suffix, "_mann_whitney_table.html"))
gtsave(out_gt, out_html)
cat(sprintf("  Tabela salva em: %s\n", out_html))

cat("\nConcluido.\n")
