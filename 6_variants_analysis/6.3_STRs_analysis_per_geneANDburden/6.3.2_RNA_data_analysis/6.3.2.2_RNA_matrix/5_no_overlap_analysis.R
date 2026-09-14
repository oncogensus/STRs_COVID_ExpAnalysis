#!/usr/bin/env Rscript
# 5_no_overlap_analysis.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Identificacao de variantes STR sem sobreposição entre case e control,
#   por intervencao (GSE). Analisa loci com outliers DBSCAN:
#     - mean_allele: (allele1 + allele2) / 2
#     - allele2_est
#   Logica: min(case) > max(control) OU min(control) > max(case)
#   Requisito: N >= 3 por grupo.
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog    intervention_outliers.tsv (saida do step 1)
#   --intervention   GSE (ex: GSE157103) ou ALL (default: ALL)
#   --out-dir        Diretorio de saida
#
# SAIDAS
#   {intervention}_no_overlap_mean_allele.csv
#   {intervention}_no_overlap_allele2.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
cat("--- No-Overlap Analysis: Case vs Control (Outlier Loci) ---\n")

df <- fread(path_str_catalog, header = TRUE, sep = "\t")
cat(sprintf("  intervention_outliers.tsv: %d linhas\n", nrow(df)))

# ==========================================
# 2. Filter by DBSCAN quality (stricter)
# ==========================================
df_filtered <- df[
  n_clusters_dbscan_global >= 1 &
  noise_ratio_dbscan_global < 0.10 &
  n_outliers_dbscan_global >= 1
]
cat(sprintf("  Apos filtro DBSCAN (clusters >= 1, noise < 10%%, outliers >= 1): %d linhas\n",

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
# Helper: identify no-overlap loci for a given metric
# ==========================================
run_no_overlap <- function(data, metric_col, min_n_per_group, label) {
  cat(sprintf("\n>>> %s (N per group >= %d) <<<\n", label, min_n_per_group))

  dt <- copy(data)
  dt[, target_metric := get(metric_col)]
  dt <- dt[!is.na(target_metric)]

  cat(sprintf("  Observacoes validas: %d\n", nrow(dt)))

  # Compute per-locus per-group statistics
  locus_stats <- dt[, .(
    median_val = median(target_metric),
    min_val = min(target_metric),
    max_val = max(target_metric),
    n_samples = .N
  ), by = .(STRs_ID, group)]

  # Pivot wider
  locus_wide <- dcast(locus_stats, STRs_ID ~ group,
                      value.var = c("median_val", "min_val", "max_val", "n_samples"),
                      fill = NA)

  # Check both groups exist
  has_case   <- "min_val_case"   in names(locus_wide)
  has_control <- "min_val_control" in names(locus_wide)

  if (!has_case || !has_control) {
    cat("  Nao ha dados para ambos os grupos.\n")
    return(data.table())
  }

  # Filter minimum sample size
  locus_wide <- locus_wide[n_samples_case >= min_n_per_group &
                           n_samples_control >= min_n_per_group]

  cat(sprintf("  Loci elegiveis (N >= %d por grupo): %d\n", min_n_per_group, nrow(locus_wide)))

  if (nrow(locus_wide) == 0) {
    cat("  Nenhum locus elegivel.\n")
    return(data.table())
  }

  # Identify no-overlap
  locus_wide[, case_higher := min_val_case > max_val_control]
  locus_wide[, control_higher := min_val_control > max_val_case]
  locus_wide[, overlap_type := case_when(
    case_higher ~ "case_higher",
    control_higher ~ "control_higher",
    TRUE ~ "overlap"
  )]

  no_overlap_ids <- locus_wide[overlap_type != "overlap"]$STRs_ID

  cat(sprintf("  Loci sem sobreposicao: %d\n", length(no_overlap_ids)))

  if (length(no_overlap_ids) == 0) {
    cat("  Nenhum locus sem sobreposicao.\n")
    return(data.table())
  }

  # Extract full observations
  final <- dt[STRs_ID %in% no_overlap_ids]
  final <- final[order(gene_name, STRs_ID)]

  cat(sprintf("  Observacoes exportadas: %d\n", nrow(final)))
  return(final)
}

# ==========================================
# 5. Run analyses
# ==========================================
# (a) Mean Allele
df_filtered[, mean_allele := (allele1_est + allele2_est) / 2]
res_mean <- run_no_overlap(df_filtered, "mean_allele", min_n_per_group = 3,
                           label = "No-Overlap: Mean Allele")

# (b) Allele 2
res_allele2 <- run_no_overlap(df_filtered, "allele2_est", min_n_per_group = 3,
                              label = "No-Overlap: Allele 2")

# ==========================================
# 6. Export results
# ==========================================
suffix <- if (intervention == "ALL") "ALL" else intervention

if (nrow(res_mean) > 0) {
  out_mean <- file.path(out_dir, paste0(suffix, "_no_overlap_mean_allele.csv"))
  write_csv2(as.data.frame(res_mean), out_mean)
  cat(sprintf("\n  Mean Allele no-overlap salvo em: %s\n", out_mean))
} else {
  cat("\n  Nenhum resultado para Mean Allele.\n")
}

if (nrow(res_allele2) > 0) {
  out_allele2 <- file.path(out_dir, paste0(suffix, "_no_overlap_allele2.csv"))
  write_csv2(as.data.frame(res_allele2), out_allele2)
  cat(sprintf("  Allele 2 no-overlap salvo em: %s\n", out_allele2))
} else {
  cat("  Nenhum resultado para Allele 2.\n")
}

cat("\nConcluido.\n")
