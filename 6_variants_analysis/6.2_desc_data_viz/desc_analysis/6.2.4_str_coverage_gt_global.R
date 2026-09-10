#!/usr/bin/env Rscript
# 6.2.4_str_coverage_gt_global.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Produces a publication-ready GT table (high-impact journal style) with the
#   GLOBAL summary of genomic STR coverage per patient (post-QC). Reads the
#   coverage_summary.tsv written by 6.2.4_str_coverage_per_patient.R.
#
# INPUTS
#   --summary   coverage_summary.tsv (GLOBAL row + per-group rows)
#   --out-dir   Output directory (default: same folder as --summary)
#
# OUTPUTS
#   table_coverage_global.html  (GT table, horizontal layout)
#   table_coverage_global.csv   (raw long-form data)
# ---------------------------------------------------------------------------
suppressMessages({
  library(data.table)
  library(gt)
})

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

args <- commandArgs(trailingOnly = TRUE)
summary_file <- get_opt(args, "--summary", NULL)
if (is.null(summary_file) || !file.exists(summary_file))
  stop("Usage: Rscript 6.2.4_str_coverage_gt_global.R --summary coverage_summary.tsv [--out-dir <dir>]")

out_dir <- get_opt(args, "--out-dir", dirname(summary_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

s <- fread(summary_file, header = TRUE, sep = "\t")
g <- s[group == "GLOBAL"]
if (nrow(g) != 1) stop("GLOBAL row not found in ", summary_file)

# ---------------------------------------------------------------------------
# FMT helpers
# ---------------------------------------------------------------------------
fmt_int   <- function(x) formatC(as.integer(round(x)), big.mark = ",", format = "d")
fmt_bp    <- function(x) formatC(as.numeric(x), digits = 0, big.mark = ",", format = "f")
fmt_mean_sd <- function(m, sd, type) {
  switch(type,
    integer = paste(fmt_int(m), "\u00b1", fmt_int(sd)),
    bp      = paste(fmt_bp(m), "\u00b1", fmt_bp(sd)),
    mb      = paste(fmt_mb(m), "\u00b1", sprintf("%.1f", as.numeric(sd) / 1e6)),
    pct     = paste(sprintf("%.4f", as.numeric(m) * 100), "\u00b1",
                    sprintf("%.4f", as.numeric(sd) * 100))
  )
}
fmt_range <- function(mn, mx, type) {
  switch(type,
    integer = paste(fmt_int(mn), "\u2013", fmt_int(mx)),
    bp      = paste(fmt_bp(mn), "\u2013", fmt_bp(mx)),
    mb      = paste(fmt_mb(mn), "\u2013", fmt_mb(mx)),
    pct     = paste(sprintf("%.4f", as.numeric(mn) * 100), "\u2013",
                    sprintf("%.4f", as.numeric(mx) * 100))
  )
}

# Internal column names are ASCII (R does not allow \u escapes in backtick names).
# Display labels with special glyphs are defined only in cols_label().
tab <- data.table(
  Variable = c(
    "STR loci per patient (n)",
    "Genomic region covered (bp)",
    "Genome coverage (%)"
  ),
  Median = c(
    fmt_int(g$n_strs_median),
    fmt_bp(g$bp_median),
    sprintf("%.4f", g$genome_frac_median * 100)
  ),
  mean_sd = c(
    fmt_mean_sd(g$n_strs_mean, g$n_strs_sd, "integer"),
    fmt_mean_sd(g$bp_mean, g$bp_sd, "bp"),
    fmt_mean_sd(g$genome_frac_mean, g$genome_frac_sd, "pct")
  ),
  Range = c(
    fmt_range(g$n_strs_min, g$n_strs_max, "integer"),
    fmt_range(g$bp_min, g$bp_max, "bp"),
    fmt_range(g$genome_frac_min, g$genome_frac_max, "pct")
  )
)

mean_sd_lab <- "Mean \u00b1 SD"

gt_tbl <- tab %>%
  gt() %>%
  cols_label(
    Variable = "Variable",
    Median = "Median",
    mean_sd = mean_sd_lab,
    Range = "Range"
  ) %>%
  tab_header(
    title = "Genomic coverage of STR loci",
    subtitle = paste0("Global summary across ", fmt_int(g$n_patients),
                      " patients (post-STRling QC)")
  ) %>%
  tab_spanner(
    label = "Per patient",
    columns = c(Median, mean_sd, Range)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = "small"),
    locations = cells_body(columns = Variable)
  ) %>%
  tab_style(
    style = cell_fill(color = "#F5F5F5"),
    locations = cells_body(columns = Variable)
  ) %>%
  tab_options(
    table.width = pct(100),
    table.border.top.style = "solid", table.border.top.width = px(2),
    table.border.bottom.style = "solid", table.border.bottom.width = px(2),
    heading.border.bottom.style = "solid", heading.border.bottom.width = px(1),
    column_labels.border.top.style = "solid", column_labels.border.top.width = px(2),
    column_labels.border.bottom.style = "solid", column_labels.border.bottom.width = px(1),
    table_body.border.bottom.style = "solid", table_body.border.bottom.width = px(1),
    data_row.padding = px(4),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(12)
  )

gtsave(gt_tbl, file.path(out_dir, "table_coverage_global.html"))
fwrite(tab, file.path(out_dir, "table_coverage_global.csv"), sep = "\t")

cat("GT table saved to:", file.path(out_dir, "table_coverage_global.html"), "\n")
cat("Raw data saved to :", file.path(out_dir, "table_coverage_global.csv"), "\n")