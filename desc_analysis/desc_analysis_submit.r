library(readr)
library(dplyr)


# Import sample-group mapping
sample_groups <- read_csv("grupos.csv")


# Remove .bam extension from sample names to match TSV format
sample_groups <- sample_groups %>%
  mutate(sample = gsub("\\.bam$", "", sample))


# Imports of STR data
df_strs <- read_delim(
  "../samples/gp_global/merged/STRs_annotated_complete.tsv",
  delim = "\t"
)


# Generate unique locus identifier
df_strs <- df_strs %>%
  mutate(
    locus = paste(chrom, start, repeat_unit, end, sep = "_")
  )


# Add group information by joining with sample_groups
df_strs <- df_strs %>%
  left_join(sample_groups, by = "sample")


# Check for unmapped samples
unmapped <- df_strs %>%
  filter(is.na(group)) %>%
  distinct(sample)


if(nrow(unmapped) > 0) {
  warning("Samples without group assignment found:")
  print(unmapped)
} else {
  cat("\nAll samples successfully mapped to groups!\n")
}


# Summary of group distribution
cat("\nGroup distribution:\n")
df_strs %>%
  count(group) %>%
  print()


# Selection of columns of interest
df_summary <- df_strs %>%
  select(region, gene_id, priority, gene_name, gene_chrom, gene_start, gene_end, 
         annotation, sample, group, locus,
         chrom, start, end, repeat_unit, allele1_est, allele2_est, depth)


# Absolute descriptive summaries
strs_by_locus <- df_summary %>% count(locus, sort = TRUE)
strs_by_region <- df_summary %>% count(region, sort = TRUE)
strs_by_repeatunit <- df_summary %>% count(repeat_unit, sort = TRUE)
strs_by_gene <- df_summary %>% count(gene_name, sort = TRUE)


# Descriptive summaries by group
strs_by_locus_group <- df_summary %>% count(locus, group, sort = TRUE)
strs_by_region_group <- df_summary %>% count(region, group, sort = TRUE)
strs_by_repeatunit_group <- df_summary %>% count(repeat_unit, group, sort = TRUE)
strs_by_gene_group <- df_summary %>% count(gene_name, group, sort = TRUE)


# Calculate allele statistics by region and group (WITH MIN AND MAX)
allele_stats_by_region_group <- df_summary %>%
  group_by(region, group) %>%
  summarise(
    mean_allele1 = mean(allele1_est, na.rm = TRUE),
    median_allele1 = median(allele1_est, na.rm = TRUE),
    sd_allele1 = sd(allele1_est, na.rm = TRUE),
    min_allele1 = min(allele1_est, na.rm = TRUE),
    max_allele1 = max(allele1_est, na.rm = TRUE),
    mean_allele2 = mean(allele2_est, na.rm = TRUE),
    median_allele2 = median(allele2_est, na.rm = TRUE),
    sd_allele2 = sd(allele2_est, na.rm = TRUE),
    min_allele2 = min(allele2_est, na.rm = TRUE),
    max_allele2 = max(allele2_est, na.rm = TRUE),
    .groups = "drop"
  )


# Create locus-to-region mapping (one region per locus)
locus_region_map <- df_summary %>%
  distinct(locus, region) %>%
  group_by(locus) %>%
  slice(1) %>%  # Keep first region if locus has multiple
  ungroup()


# Create repeat_unit summary with groups info (NOT separated by group)
strs_by_repeatunit_combo <- strs_by_repeatunit_group %>%
  group_by(repeat_unit) %>%
  summarise(
    n_case = sum(n[group == "case"], na.rm = TRUE),
    n_control = sum(n[group == "control"], na.rm = TRUE),
    total = sum(n),
    groups_present = paste(sort(unique(group)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(total))


# strs_by_locus_combo: DETAILED FORMAT - one row per variant
strs_by_locus_combo <- df_summary %>%
  select(locus, sample, group, region, chrom, start, end, repeat_unit, 
         allele1_est, allele2_est, depth) %>%
  arrange(locus, group, sample)


strs_by_region_combo <- strs_by_region_group %>%
  left_join(strs_by_region %>% rename(total = n), by = "region") %>%
  left_join(allele_stats_by_region_group, by = c("region", "group"))


strs_by_gene_combo <- left_join(strs_by_gene_group,
                                 strs_by_gene %>% rename(total = n),
                                 by = "gene_name")


# Define output directory
out_dir <- "results/"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Save all combos as CSV
write_csv(strs_by_locus_combo, file.path(out_dir, "strs_by_locus_combo.csv"))
write_csv(strs_by_region_combo, file.path(out_dir, "strs_by_region_combo.csv"))
write_csv(strs_by_repeatunit_combo, file.path(out_dir, "strs_by_repeatunit_combo.csv"))
write_csv(strs_by_gene_combo, file.path(out_dir, "strs_by_gene_combo.csv"))


# Calculate mean allele and extension size
df_strs <- df_strs %>%
  mutate(
    mean_allele = rowMeans(cbind(allele1_est, allele2_est), na.rm = TRUE),
    extension_size = nchar(repeat_unit) * mean_allele
  )


# Descriptive statistics by group
stats_group_repeat <- df_strs %>%
  group_by(group) %>%
  summarise(
    Mean_Allele1   = mean(allele1_est, na.rm = TRUE),
    SD_Allele1     = sd(allele1_est, na.rm = TRUE),
    Q1_Allele1     = quantile(allele1_est, 0.25, na.rm = TRUE),
    Median_Allele1 = median(allele1_est, na.rm = TRUE),
    Q3_Allele1     = quantile(allele1_est, 0.75, na.rm = TRUE),
    
    Mean_Allele2   = mean(allele2_est, na.rm = TRUE),
    SD_Allele2     = sd(allele2_est, na.rm = TRUE),
    Q1_Allele2     = quantile(allele2_est, 0.25, na.rm = TRUE),
    Median_Allele2 = median(allele2_est, na.rm = TRUE),
    Q3_Allele2     = quantile(allele2_est, 0.75, na.rm = TRUE),
    
    Mean_ExtensionSize   = mean(extension_size, na.rm = TRUE),
    SD_ExtensionSize     = sd(extension_size, na.rm = TRUE),
    Q1_ExtensionSize     = quantile(extension_size, 0.25, na.rm = TRUE),
    Median_ExtensionSize = median(extension_size, na.rm = TRUE),
    Q3_ExtensionSize     = quantile(extension_size, 0.75, na.rm = TRUE),
    
    Number_Observations = n(),
    .groups = "drop"
  )


# Export results
write_csv(stats_group_repeat, file.path(out_dir, "stats_group_repeat_alleles_meanallele_extension_quartiles_sd.csv"))


cat("\nAnalysis completed successfully!\n")
cat("Results saved in:", out_dir, "\n")
