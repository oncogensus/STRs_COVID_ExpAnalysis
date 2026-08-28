library(data.table)
library(dplyr)
library(stringr)

# ==========================================
# STEP 1: Load Integrated Dataset
# ==========================================
cat("\n", strrep("=", 42), "\n")
cat("STEP 1: Loading Data\n")
cat(strrep("=", 42), "\n")

dataset_path <- "../../../samples/STRs_analysis_dataset.tsv"
df <- fread(dataset_path)

# Creating a unique locus ID (chrom_start_unit_end)
df <- df %>%
  mutate(locus = paste(chrom, start, repeat_unit, end, sep = "_"))

# Quick Group Check
cat("\nGroup distribution:\n")
print(df %>% count(group))

# ==========================================
# STEP 2: Calculate Derived Metrics
# ==========================================
df <- df %>%
  mutate(
    mean_allele = (allele1_est + allele2_est) / 2,
    extension_size = nchar(repeat_unit) * mean_allele
  )

# ==========================================
# STEP 3: Generate Detailed Summaries
# ==========================================
cat("\nGenerating Detailed Summaries...\n")

# 3A. By Locus (Detailed variant-level view)
strs_by_locus_combo <- df %>%
  select(locus, sample_id, group, region, gene_name, chrom, start, end, 
         repeat_unit, allele1_est, allele2_est, depth) %>%
  arrange(locus, group, sample_id)

# 3B. By Region (Including stats)
strs_by_region_combo <- df %>%
  group_by(region, group) %>%
  summarise(
    n_per_group  = n(),
    mean_allele1 = mean(allele1_est, na.rm = TRUE),
    sd_allele1   = sd(allele1_est, na.rm = TRUE),
    min_allele1  = min(allele1_est, na.rm = TRUE),
    max_allele1  = max(allele1_est, na.rm = TRUE),
    mean_allele2 = mean(allele2_est, na.rm = TRUE),
    sd_allele2   = sd(allele2_est, na.rm = TRUE),
    min_allele2  = min(allele2_est, na.rm = TRUE),
    max_allele2  = max(allele2_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(df %>% count(region, name = "total_region"), by = "region")

# 3C. By Repeat Unit (With Case/Control counts)
strs_by_repeatunit_combo <- df %>%
  group_by(repeat_unit) %>%
  summarise(
    n_case         = sum(group == "case", na.rm = TRUE),
    n_control      = sum(group == "control", na.rm = TRUE),
    total          = n(),
    groups_present = paste(sort(unique(group)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(total))

# 3D. By Gene (New addition from previous script)
strs_by_gene_combo <- df %>%
  group_by(gene_name, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(df %>% count(gene_name, name = "total_gene"), by = "gene_name") %>%
  arrange(desc(total_gene))

# ==========================================
# STEP 4: Publication Stats (GT-Ready)
# ==========================================
cat("Calculating Publication Stats...\n")

stats_group_repeat <- df %>%
  group_by(group) %>%
  summarise(
    across(c(allele1_est, allele2_est, extension_size), 
           list(
             Mean   = ~mean(.x, na.rm = TRUE),
             SD     = ~sd(.x, na.rm = TRUE),
             Median = ~median(.x, na.rm = TRUE),
             Q1     = ~quantile(.x, 0.25, na.rm = TRUE),
             Q3     = ~quantile(.x, 0.75, na.rm = TRUE)
           ),
           .names = "{.fn}_{.col}"),
    Number_Observations = n(),
    .groups = "drop"
  )

# ==========================================
# STEP 5: Exporting Results
# ==========================================
out_dir <- "results/"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fwrite(strs_by_locus_combo,      file.path(out_dir, "strs_by_locus_combo.csv"))
fwrite(strs_by_region_combo,     file.path(out_dir, "strs_by_region_combo.csv"))
fwrite(strs_by_repeatunit_combo, file.path(out_dir, "strs_by_repeatunit_combo.csv"))
fwrite(strs_by_gene_combo,       file.path(out_dir, "strs_by_gene_combo.csv"))
fwrite(stats_group_repeat,       file.path(out_dir, "stats_group_repeat_full_metrics.csv"))

cat("\nAnalysis completed successfully! Results in:", out_dir, "\n")