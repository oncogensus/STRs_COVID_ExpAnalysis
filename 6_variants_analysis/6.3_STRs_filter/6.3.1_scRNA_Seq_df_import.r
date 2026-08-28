library(dplyr)
library(tidyr)
library(readr)

# =======================================================
# 1. LOAD DATASETS
# =======================================================
cat("--- Loading Datasets ---\n")

# STR Dataset (Your DBSCAN results)
str_path <- "../../samples/STRs_analysis_dataset.tsv"
df_strs <- read_tsv(str_path)

# =======================================================
# 2. FILTER STR OUTLIERS (Relaxed Quality)
# =======================================================
cat("--- Filtering STR Outliers (Noise <= 0.10) ---\n")

df_str_candidates <- df_strs %>%
  filter(n_clusters > 0) %>%
  filter(noise_ratio <= 0.10) %>%
  filter(n_outliers >= 1) %>%
  mutate(
    outlier_samples = as.character(outlier_samples),
    outlier_residuals = as.character(outlier_residuals)
  ) %>%
  separate_rows(
    outlier_samples, 
    outlier_residuals, 
    sep = ";"
  ) %>%
  mutate(abs_res = abs(as.numeric(outlier_residuals))) %>%
  group_by(STRs_ID) %>%
  slice_max(abs_res, n = 1, with_ties = FALSE) %>%
  ungroup()

# =======================================================
# 3. PROCESS EACH GSE FOLDER SEPARATELY
# =======================================================
cat("--- Crossing STR Outliers with Single-Cell DEGs per GSE ---\n")

sc_base <- "../../6_scovid_data"
gse_folders <- list.dirs(sc_base, recursive = FALSE)

if (!dir.exists("results")) dir.create("results")

gse_to_tissue <- c(
  "GSE157344" = "lung",
  "GSE159812" = "brain"
)

all_results <- list()

for (gse in gse_folders) {
  gse_name <- basename(gse)
  cat(sprintf("\n=== Processing %s ===\n", gse_name))

  csv_files <- list.files(gse, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("_total\\.csv$", csv_files)]
  cat(sprintf("  Found %d CSV files (%d after excluding totals)\n", length(csv_files), length(csv_files)))

  df_gse <- csv_files %>%
    lapply(function(f) {
      df <- read_csv(f, show_col_types = FALSE)
      cell_type <- f %>%
        basename() %>%
        tools::file_path_sans_ext() %>%
        stringr::str_remove(paste0("^", gse_name, "_"))
      df$cell_type <- cell_type
      df
    }) %>%
    bind_rows()

  overlap_analysis <- df_str_candidates %>%
    inner_join(df_gse, by = c("gene_name" = "Gene Symbol")) %>%
    mutate(GEO = gse_name, source_tissue = gse_to_tissue[gse_name]) %>%
    select(
      gene_name,
      cell_type,
      LogFC,
      Pvalue,
      abs_res,
      STRs_ID,
      outlier_samples,
      noise_ratio,
      outlier_residuals,
      n_clusters,
      depth,
      allele1_est,
      allele2_est,
      GEO,
      source_tissue,
      group
    ) %>%
    distinct() %>%
    arrange(desc(abs_res))

  out_file <- sprintf("results/STR_vs_scRNA_overlap_%s.csv", gse_name)
  write_csv(overlap_analysis, out_file)

  all_results[[gse_name]] <- overlap_analysis

  cat(sprintf("  Genes found: %d\n", n_distinct(overlap_analysis$gene_name)))
  cat(sprintf("  Interactions: %d\n", nrow(overlap_analysis)))
  cat(sprintf("  Saved: %s\n", out_file))
}

# =======================================================
# 4. UNIFIED OUTPUT WITH SOURCE TISSUE
# =======================================================
cat("\n--- Generating unified output ---\n")

unified <- bind_rows(all_results)
write_csv(unified, "results/STR_vs_scRNA_overlap_unified.csv")
cat(sprintf("Unified: results/STR_vs_scRNA_overlap_unified.csv (%d rows)\n", nrow(unified)))

cat("\n=== All done ===\n")
