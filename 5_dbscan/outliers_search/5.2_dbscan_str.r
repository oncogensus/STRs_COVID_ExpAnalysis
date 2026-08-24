#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(dbscan)

# ---------------------------
# Configuration
# ---------------------------
minpts <- 2                     # base value for minPts calculation
eps_default <- 2                 # not directly used; kept for reference
input_file <- "../norm_test/STRs_normalized_residuals.tsv"
output_file <- "results_dbscan/outliers_per_str.tsv"

# ---------------------------
# Read consolidated data
# ---------------------------
cat("Reading file:", input_file, "\n")
dat <- fread(input_file, header = TRUE, sep = "\t", data.table = FALSE)

# Check for required columns
required_cols <- c("STRs_ID", "allele2_residuals", "sample_id", "group")
if (!all(required_cols %in% colnames(dat))) {
  stop("Input file must contain the following columns: ", paste(required_cols, collapse = ", "))
}

# ---------------------------
# Prepare output structure
# ---------------------------
# Get unique STR identifiers
str_list <- unique(dat$STRs_ID)
cat("Total STRs:", length(str_list), "\n")

# Data frame to store results (one row per STR)
results <- data.frame(
  STRs_ID = str_list,
  n_samples = NA,
  n_samples_valid = NA,  # NOVA: amostras válidas (sem NA)
  n_outliers = NA,
  outlier_samples = NA,
  outlier_residuals = NA,
  n_clusters = NA,
  noise_ratio = NA,
  stringsAsFactors = FALSE
)

# ---------------------------
# Process each STR individually
# ---------------------------
for (i in seq_along(str_list)) {
  if (i %% 100 == 0) cat("Processing STR", i, "of", length(str_list), "\n")
  
  str_id <- str_list[i]
  
  # Subset data for this STR
  dat_str <- subset(dat, STRs_ID == str_id)
  
  # Number of samples for this STR (total, including NAs)
  n_total <- nrow(dat_str)
  results[i, "n_samples"] <- n_total
  
  # Remove rows with NA residuals
  dat_str <- dat_str[!is.na(dat_str$allele2_residuals), ]
  n_valid <- nrow(dat_str)
  results[i, "n_samples_valid"] <- n_valid
  
  if (n_valid < 3) {
    next  # Skip STRs with insufficient valid samples
  }
  
  # Extract residuals and sample IDs
  residuals <- dat_str$allele2_residuals
  sample_ids <- dat_str$sample_id
  
  # ---------------------------
  # Calculate DBSCAN parameters
  # ---------------------------
  # minPts: based on valid sample size
  minPts <- ceiling(log2(minpts * n_valid))
  minPts <- max(2, minPts)      # ensure at least 2
  
  # eps: based on residual spread
  eps_value <- as.numeric(diff(quantile(residuals, c(0.05, 0.95), na.rm = TRUE)))
  eps_value <- max(eps_value, 1e-6)   # avoid zero
  
  # ---------------------------
  # Run DBSCAN
  # ---------------------------
  X <- matrix(residuals, ncol = 1)
  scan <- dbscan(X, eps = eps_value, minPts = minPts)
  
  # ---------------------------
  # Determine cutoff and outliers
  # ---------------------------
  clusters <- scan$cluster
  unique_clusters <- unique(clusters)
  
  if (length(unique_clusters) == 1 || sum(clusters == 0) == 0) {
    # Only one cluster or no noise → no outliers (cutoff = Inf)
    cutoff <- Inf
  } else {
    # Largest residual among non‑noise points (cluster != 0)
    non_noise_residuals <- residuals[clusters != 0]
    cutoff <- max(non_noise_residuals, na.rm = TRUE)
    cutoff <- ifelse(cutoff < 2, 2, cutoff)   # original minimum of 2
  }
  
  # Identify outliers (residuals above cutoff)
  outlier_idx <- which(residuals > cutoff)
  n_outliers <- length(outlier_idx)
  
  # ---------------------------
  # Store results for this STR
  # ---------------------------
  results[i, "n_outliers"] <- n_outliers
  results[i, "outlier_samples"] <- paste(sample_ids[outlier_idx], collapse = ";")
  results[i, "outlier_residuals"] <- paste(residuals[outlier_idx], collapse = ";")
  results[i, "n_clusters"] <- length(unique_clusters[unique_clusters != 0])
  results[i, "noise_ratio"] <- sum(clusters == 0) / n_valid
}

# ---------------------------
# Save results
# ---------------------------
fwrite(results, output_file, sep = "\t", row.names = FALSE)
cat("Results saved to:", output_file, "\n")