#!/usr/bin/env Rscript
# 7_merged_table_allele_stats.R
# ---------------------------------------------------------------------------
# Generates a single transposed (wide) table with descriptive statistics of
# allele2_est per genomic region + overall, comparing Case vs Control.
#
# Inputs:
#   results/strs_by_locus_combo.csv  (raw data: locus, sample_id, group,
#                                     region, allele2_est, ...)
#
# Outputs:
#   results/table_allele_stats_merged.html
#   results/table_allele_stats_merged.csv
#
# Usage:
#   Rscript 7_merged_table_allele_stats.R
#   or as a cell in 2_desc_analysis.ipynb
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(gt)
})

# ==========================================
# 1. Load raw data
# ==========================================
df <- read_csv("results/strs_by_locus_combo.csv", show_col_types = FALSE)

cat(sprintf("Data loaded: %d rows, %d unique regions\n",
            nrow(df), n_distinct(df$region)))

# ==========================================
# 2. Compute statistics per region + group
# ==========================================
region_stats <- df %>%
  group_by(region, group) %>%
  summarise(
    N        = n(),
    Mean     = mean(allele2_est, na.rm = TRUE),
    SD       = sd(allele2_est, na.rm = TRUE),
    Median   = median(allele2_est, na.rm = TRUE),
    Q1       = quantile(allele2_est, 0.25, na.rm = TRUE),
    Q3       = quantile(allele2_est, 0.75, na.rm = TRUE),
    Min      = min(allele2_est, na.rm = TRUE),
    Max      = max(allele2_est, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  # Readable region names
  mutate(
    region = case_when(
      region == "CDS"              ~ "CDS",
      region == "five_prime_utr"   ~ "5' UTR",
      region == "three_prime_utr"  ~ "3' UTR",
      region == "intergenic"       ~ "Intergenic",
      region == "intron"           ~ "Intron",
      region == "promoter"         ~ "Promoter",
      region == "non_coding_exons" ~ "Non-coding Exons",
      region == "others"           ~ "Others",
      TRUE ~ region
    )
  )

# ==========================================
# 3. Global statistics (Overall)
# ==========================================
overall_stats <- df %>%
  group_by(group) %>%
  summarise(
    N        = n(),
    Mean     = mean(allele2_est, na.rm = TRUE),
    SD       = sd(allele2_est, na.rm = TRUE),
    Median   = median(allele2_est, na.rm = TRUE),
    Q1       = quantile(allele2_est, 0.25, na.rm = TRUE),
    Q3       = quantile(allele2_est, 0.75, na.rm = TRUE),
    Min      = min(allele2_est, na.rm = TRUE),
    Max      = max(allele2_est, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(region = "Overall")

# ==========================================
# 4. Join and transpose to wide (Option B)
# ==========================================
all_stats <- bind_rows(overall_stats, region_stats)

# Transpose: metrics become columns with _Case / _Control suffix
wide_case <- all_stats %>%
  filter(group == "case") %>%
  select(-group) %>%
  rename_with(~ paste0(.x, "_Case"), -c(region))

wide_control <- all_stats %>%
  filter(group == "control") %>%
  select(-group) %>%
  rename_with(~ paste0(.x, "_Control"), -c(region))

merged <- full_join(wide_case, wide_control, by = "region") %>%
  # Desired column order
  select(
    region,
    N_Case, N_Control,
    Mean_Case, Mean_Control,
    SD_Case, SD_Control,
    Median_Case, Median_Control,
    Q1_Case, Q1_Control,
    Q3_Case, Q3_Control,
    Min_Case, Min_Control,
    Max_Case, Max_Control
  ) %>%
  # Order regions (Overall first)
  mutate(region = factor(region, levels = c(
    "Overall", "Promoter", "5' UTR", "CDS", "3' UTR",
    "Intron", "Non-coding Exons", "Intergenic", "Others"
  ))) %>%
  arrange(region)

cat(sprintf("Final table: %d rows x %d columns\n", nrow(merged), ncol(merged)))

# ==========================================
# 5. Generate GT table
# ==========================================
gt_table <- merged %>%
  gt(rowname_col = "region") %>%

  # Spanners for Case and Control
  tab_spanner(
    label = html("<b>Case</b>"),
    columns = ends_with("_Case")
  ) %>%
  tab_spanner(
    label = html("<b>Control</b>"),
    columns = ends_with("_Control")
  ) %>%

  # Column labels (remove suffix)
  cols_label(
    N_Case       = "N",
    N_Control    = "N",
    Mean_Case    = "Mean",
    Mean_Control = "Mean",
    SD_Case      = "SD",
    SD_Control   = "SD",
    Median_Case  = "Median",
    Median_Control = "Median",
    Q1_Case      = "Q1",
    Q1_Control   = "Q1",
    Q3_Case      = "Q3",
    Q3_Control   = "Q3",
    Min_Case     = "Min",
    Min_Control  = "Min",
    Max_Case     = "Max",
    Max_Control  = "Max"
  ) %>%

  # Number formatting
  fmt_number(
    columns = where(is.numeric) & !ends_with("N_Case") & !ends_with("N_Control"),
    decimals = 2,
    use_seps = TRUE
  ) %>%
  fmt_number(
    columns = c(N_Case, N_Control),
    decimals = 0,
    use_seps = TRUE
  ) %>%

  # Title
  tab_header(
    title = "Descriptive Statistics of Major Allele Size (allele2_est)",
    subtitle = "By Genomic Region — Case vs Control"
  ) %>%

  # Footnote
  tab_footnote(
    footnote = "SD: Standard Deviation; Q1: First Quartile; Q3: Third Quartile",
    locations = cells_column_labels(columns = c(SD_Case, SD_Control))
  ) %>%

  # Style
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_fill(color = "grey95"),
    locations = cells_row_groups()
  ) %>%

  # Borders
  tab_options(
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    heading.border.bottom.style = "solid",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table_body.border.bottom.style = "solid",
    row_group.border.top.style = "solid",
    row_group.border.bottom.style = "solid",
    table.width = pct(100)
  )

# ==========================================
# 6. Save
# ==========================================
gtsave(gt_table, "results/table_allele_stats_merged.html")
write_csv(merged, "results/table_allele_stats_merged.csv")

cat("\n=== Outputs ===\n")
cat("  results/table_allele_stats_merged.html\n")
cat("  results/table_allele_stats_merged.csv\n")
cat("=== END ===\n")