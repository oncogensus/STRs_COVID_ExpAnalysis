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

# scRNA-seq DEGs (The COVID-19 single-cell data)
sc_path <- "../../6_scovid_data/smerged_covid.csv"
df_scRNA <- read_csv(sc_path)

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
  # CORREÇÃO 1: Se GEO, sample_id_clean e outros metadados vieram 
  # concatenados com ";" no df_strs, eles DEVEM ser incluídos aqui!
  separate_rows(
    outlier_samples, 
    outlier_residuals, 
    # DESCOMENTE as linhas abaixo se elas também tiverem ";" no seu dataset original:
    # GEO, 
    # sample_id_clean, 
    # source_tissue, 
    # group,
    sep = ";"
  ) %>%
  mutate(abs_res = abs(as.numeric(outlier_residuals))) %>%
  # Keep the best outlier per locus
  group_by(STRs_ID) %>%
  slice_max(abs_res, n = 1, with_ties = FALSE) %>%
  ungroup()

# =======================================================
# 3. INTERSECT STR GENES WITH DEGs
# =======================================================
# =======================================================
# 3. INTERSECT STR GENES WITH DEGs
# =======================================================
cat("--- Crossing STR Outliers with Single-Cell DEGs ---\n")

overlap_analysis <- df_str_candidates %>%
  inner_join(df_scRNA, by = c("gene_name" = "Gene Symbol")) %>%
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
    GEO = GEO_ID,          # <- CORREÇÃO: Mapeia a coluna populada (GEO_ID) para o nome GEO
    source_tissue,
    group
  ) %>%
  distinct() %>%           # Garante que não haverá linhas idênticas duplicadas
  arrange(desc(abs_res))

# =======================================================
# 4. RESULTS & EXPORT
# =======================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total STR-affected genes found in scRNA-seq data:", n_distinct(overlap_analysis$gene_name), "\n")
cat("Total cell-type specific interactions found:", nrow(overlap_analysis), "\n")

if (nrow(overlap_analysis) > 0) {
  print(head(overlap_analysis, 15))
  
  if (!dir.exists("results")) dir.create("results")
  write_csv(overlap_analysis, "results/STR_vs_scRNA_overlap.csv")
  cat("\n✓ Results saved to 'results/STR_vs_scRNA_overlap.csv'\n")
} else {
  cat("\n⚠️ No overlap found between STR outliers and DEGs.\n")
}