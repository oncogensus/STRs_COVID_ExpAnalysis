#!/usr/bin/env Rscript
# 4_mann_whitney_allele.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Mann-Whitney U test: case vs. control para STRs com outliers DBSCAN,
#   por intervencao (GSE). Dois testes:
#     (a) allele2_est com N > 3 por grupo
#     (b) mean_allele = (allele1 + allele2) / 2 com N >= 3 por grupo
#   Aplica FDR (Benjamini-Hochberg) para correcao de multiplos testes.
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
cat(sprintf("  intervention_outliers.tsv: %d linhas\n", nrow(df)))

# ==========================================
# 2. Filter by DBSCAN quality (stricter than file default)
# ==========================================
df_filtered <- df[
  n_clusters_dbscan_global > 1 &
  noise_ratio_dbscan_global < 0.10 &
  n_outliers_dbscan_global > 1
]
cat(sprintf("  Apos filtro DBSCAN (clusters > 1, noise < 10%%, outliers > 1): %d linhas\n",
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

  # Filter loci with sufficient samples per group
  eligible_loci <- dt[, .N, by = .(STRs_ID, group)]
  eligible_loci <- dcast(eligible_loci, STRs_ID ~ group, value.var = "N", fill = 0)
  eligible_loci <- eligible_loci[case >= min_n_per_group & control >= min_n_per_group]
  eligible_ids <- eligible_loci$STRs_ID

  cat(sprintf("  Loci elegiveis: %d\n", length(eligible_ids)))

  if (length(eligible_ids) == 0) {
    cat("  Nenhum locus elegivel.\n")
    return(data.table())
  }

  # Run Wilcoxon test per locus
  dt_eligible <- dt[STRs_ID %in% eligible_ids]

  results <- dt_eligible %>%
    group_by(STRs_ID) %>%
    wilcox_test(target_metric ~ group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    as.data.table()

  sig_results <- results[p.adj < 0.05]
  sig_ids <- sig_results$STRs_ID

  cat(sprintf("  Significativos (p.adj < 0.05): %d\n", length(sig_ids)))

  if (length(sig_ids) == 0) {
    cat("  Nenhum locus significativo apos FDR.\n")
    return(data.table())
  }

  # Extract full observations for significant loci
  final <- dt[STRs_ID %in% sig_ids]
  final <- merge(final,
                 sig_results[, .(STRs_ID, p, p.adj, p.adj.signif)],
                 by = "STRs_ID", all.x = TRUE)
  final <- final[order(p.adj, STRs_ID)]

  cat(sprintf("  Observacoes exportadas: %d\n", nrow(final)))
  return(final)
}

# ==========================================
# 5. Run tests
# ==========================================
# (a) Allele 2
res_allele2 <- run_mann_whitney(df_filtered, "allele2_est", min_n_per_group = 4,
                                label = "Mann-Whitney: Allele 2")

# (b) Mean Allele
df_filtered[, mean_allele := (allele1_est + allele2_est) / 2]
res_mean <- run_mann_whitney(df_filtered, "mean_allele", min_n_per_group = 3,
                             label = "Mann-Whitney: Mean Allele")

# ==========================================
# 6. Export results
# ==========================================
suffix <- if (intervention == "ALL") "ALL" else intervention

if (nrow(res_allele2) > 0) {
  out_allele2 <- file.path(out_dir, paste0(suffix, "_allele2_mw.csv"))
  write_csv2(as.data.frame(res_allele2), out_allele2)
  cat(sprintf("\n  Allele2 MW salvo em: %s\n", out_allele2))
} else {
  cat("\n  Nenhum resultado significativo para Allele 2.\n")
}

if (nrow(res_mean) > 0) {
  out_mean <- file.path(out_dir, paste0(suffix, "_mean_allele_mw.csv"))
  write_csv2(as.data.frame(res_mean), out_mean)
  cat(sprintf("  Mean Allele MW salvo em: %s\n", out_mean))
} else {
  cat("  Nenhum resultado significativo para Mean Allele.\n")
}

cat("\nConcluido.\n")
