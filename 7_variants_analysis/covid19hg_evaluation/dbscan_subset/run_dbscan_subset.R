#!/usr/bin/env Rscript
# run_dbscan_subset.R
# Re-roda o DBSCAN (mesmos parametros de 5_dbscan/outliers_search/5.2_dbscan_str.r)
# sobre o SUBSET de STRs localizados nos genes sugestivos do COVID-19 HG.
#
# Uso:
#   Rscript run_dbscan_subset.R <input_residuals.tsv> <output.tsv>

# Garante que a biblioteca de usuario (R_LIBS_USER) esteja no path e instala
# data.table/dbscan caso nao existam no ambiente (ex.: env `igv` sem os pacotes).
rlibs <- Sys.getenv("R_LIBS_USER")
if (rlibs != "") .libPaths(c(rlibs, .libPaths()))
for (pkg in c("data.table", "dbscan")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org",
                     lib = if (rlibs != "") rlibs else NULL)
  }
}
library(data.table)
library(dbscan)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript run_dbscan_subset.R <input> <output>")
input_file <- args[1]
output_file <- args[2]

minpts <- 2
dat <- fread(input_file, header = TRUE, sep = "\t", data.table = FALSE)

required_cols <- c("STRs_ID", "allele2_residuals", "sample_id", "group")
if (!all(required_cols %in% colnames(dat)))
  stop("Input deve conter: ", paste(required_cols, collapse = ", "))

str_list <- unique(dat$STRs_ID)
cat("Total STRs (subset):", length(str_list), "\n")

results <- data.frame(
  STRs_ID = str_list,
  n_samples = NA,
  n_samples_valid = NA,
  n_outliers = NA,
  outlier_samples = NA,
  outlier_residuals = NA,
  n_clusters = NA,
  noise_ratio = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(str_list)) {
  if (i %% 100 == 0) cat("Processing STR", i, "of", length(str_list), "\n")
  str_id <- str_list[i]
  dat_str <- subset(dat, STRs_ID == str_id)
  n_total <- nrow(dat_str)
  results[i, "n_samples"] <- n_total
  dat_str <- dat_str[!is.na(dat_str$allele2_residuals), ]
  n_valid <- nrow(dat_str)
  results[i, "n_samples_valid"] <- n_valid
  if (n_valid < 3) next
  residuals <- dat_str$allele2_residuals
  sample_ids <- dat_str$sample_id

  minPts <- ceiling(log2(minpts * n_valid))
  minPts <- max(2, minPts)
  eps_value <- as.numeric(diff(quantile(residuals, c(0.05, 0.95), na.rm = TRUE)))
  eps_value <- max(eps_value, 1e-6)

  X <- matrix(residuals, ncol = 1)
  scan <- dbscan(X, eps = eps_value, minPts = minPts)
  clusters <- scan$cluster
  unique_clusters <- unique(clusters)

  if (length(unique_clusters) == 1 || sum(clusters == 0) == 0) {
    cutoff <- Inf
  } else {
    non_noise_residuals <- residuals[clusters != 0]
    cutoff <- max(non_noise_residuals, na.rm = TRUE)
    cutoff <- ifelse(cutoff < 2, 2, cutoff)
  }

  outlier_idx <- which(residuals > cutoff)
  n_outliers <- length(outlier_idx)
  results[i, "n_outliers"] <- n_outliers
  results[i, "outlier_samples"] <- paste(sample_ids[outlier_idx], collapse = ";")
  results[i, "outlier_residuals"] <- paste(residuals[outlier_idx], collapse = ";")
  results[i, "n_clusters"] <- length(unique_clusters[unique_clusters != 0])
  results[i, "noise_ratio"] <- sum(clusters == 0) / n_valid
}

fwrite(results, output_file, sep = "\t", row.names = FALSE)
cat("Results saved to:", output_file, "\n")
