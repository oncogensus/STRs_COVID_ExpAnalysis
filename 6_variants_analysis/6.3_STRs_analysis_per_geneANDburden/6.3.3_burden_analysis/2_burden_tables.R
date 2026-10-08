#!/usr/bin/env Rscript
# 2_burden_tables.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Builds publication-ready supplementary tables (gt HTML) from the burden
#   analysis outputs (step 1_burden_analysis.R). Assembles Supplementary
#   Table S11 with four panels:
#     Panel A - Global relative burden (Mann-Whitney)
#     Panel B - Firth logistic regression of global relative burden
#     Panel C - Relative burden within DEGs (global + per intervention)
#     Panel D - Relative burden per genomic region
#
# INPUTS (via command-line arguments)
#   --results-dir   Directory with burden_*.csv (default: results)
#   --out-dir       Output directory for tables (default: <results-dir>/tables)
#
# OUTPUTS
#   Supplementary_Table_S11.html      All four panels in one document
#   Supplementary_Table_S11_A.html    Panel A
#   Supplementary_Table_S11_B.html    Panel B
#   Supplementary_Table_S11_C.html    Panel C
#   Supplementary_Table_S11_D.html    Panel D
#
# ENVIRONMENT
#   r_enrich_env (micromamba): data.table, gt
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(gt)
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

results_dir <- parse_arg("--results-dir", "results")
out_dir     <- parse_arg("--out-dir", file.path(results_dir, "tables"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# Helpers
# ==========================================
fmt_pval <- function(x) {
  ifelse(is.na(x), "Not estimable",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                ifelse(x < 0.05, sprintf("%.3f", x), sprintf("%.2f", x))))
}

fmt_num <- function(x, d = 2) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))

fmt_or_ci <- function(or, lo, hi, d = 2) {
  ifelse(is.na(or), "-", sprintf("%s (%s\u2013%s)", fmt_num(or, d), fmt_num(lo, d), fmt_num(hi, d)))
}

read_csv_safe <- function(f) {
  if (file.exists(f)) return(fread(f))
  cat(sprintf("  [WARN] missing: %s\n", f))
  NULL
}

# Common publication style.
style_gt <- function(g) {
  g %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_options(
      table.font.size = px(13),
      table.width = pct(100),
      heading.align = "center",
      heading.title.font.size = px(15),
      heading.subtitle.font.size = px(12),
      column_labels.font.weight = "bold",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      table_body.border.bottom.style = "solid",
      data_row.padding = px(5)
    )
}

legend_notes <- list(
  md("*OR*, odds ratio; *CI*, confidence interval; *DEG*, differentially expressed gene; *IQR*, interquartile range; *PC1*, first ancestry principal component."),
  md("Relative burden = proportion of DBSCAN-identified outlier STRs per individual. Logistic regression was fitted using Firth's penalized likelihood. All *P*-values are nominal and exploratory (no multiple-testing correction).")
)

cat("--- Supplementary Table S11 (burden analyses) ---\n")

# ==========================================
# Panel A - Global relative burden (Mann-Whitney)
# ==========================================
gtA <- NULL
gmw <- read_csv_safe(file.path(results_dir, "burden_global_mw.csv"))
if (!is.null(gmw)) {
  pa <- data.table(
    Group = c("Survivors", "Fatal COVID-19 cases"),
    N = c(gmw$n_survivors[1], gmw$n_fatal[1]),
    `Relative burden, median [IQR]` = c(gmw$survivors_median_IQR[1], gmw$fatal_median_IQR[1])
  )
  gtA <- pa %>%
    gt() %>%
    tab_header(title = md("**Supplementary Table S11 \u2014 Panel A. Global relative burden of outlier STRs**")) %>%
    cols_label(N = "N", `Relative burden, median [IQR]` = md("Relative burden, median [IQR]")) %>%
    style_gt() %>%
    tab_source_note(source_note = md(sprintf("Mann\u2013Whitney U *P* = %s.", fmt_pval(gmw$p[1])))) %>%
    tab_source_note(source_note = legend_notes[[1]])
}

# ==========================================
# Panel B - Firth logistic regression
# ==========================================
gtB <- NULL
fb <- read_csv_safe(file.path(results_dir, "burden_firth.csv"))
if (!is.null(fb)) {
  pb <- data.table(
    Predictor = fb$predictor,
    `OR (95% CI)` = fmt_or_ci(fb$OR, fb$CI_low, fb$CI_high),
    `P (nominal)` = fmt_pval(fb$p)
  )
  gtB <- pb %>%
    gt() %>%
    tab_header(title = md("**Supplementary Table S11 \u2014 Panel B. Firth logistic regression of global relative burden on fatality**")) %>%
    cols_label(Predictor = "Predictor", `OR (95% CI)` = md("OR (95% CI)"), `P (nominal)` = md("*P* (nominal)")) %>%
    style_gt() %>%
    tab_source_note(source_note = md("Outcome: fatal COVID-19 (1) vs survivor (0). Predictors: relative burden per 1 percentage point, age per year, sex (male vs female), PC1 per standard deviation.")) %>%
    tab_source_note(source_note = legend_notes[[2]])
}

# ==========================================
# Panel C - Relative burden within DEGs
# ==========================================
gtC <- NULL
dc <- read_csv_safe(file.path(results_dir, "burden_deg_mw.csv"))
if (!is.null(dc)) {
  pc <- data.table(
    Context = dc$context,
    `N (survivors/fatal)` = paste0(dc$n_survivors, "/", dc$n_fatal),
    `Survivors, median [IQR]` = dc$survivors_median_IQR,
    `Fatal, median [IQR]` = dc$fatal_median_IQR,
    `P (nominal)` = fmt_pval(dc$p)
  )
  gtC <- pc %>%
    gt() %>%
    tab_header(title = md("**Supplementary Table S11 \u2014 Panel C. Relative burden within DEGs**")) %>%
    cols_label(Context = "Context", `N (survivors/fatal)` = "N (survivors/fatal)",
               `Survivors, median [IQR]` = "Survivors, median [IQR]",
               `Fatal, median [IQR]` = "Fatal, median [IQR]",
               `P (nominal)` = md("*P* (nominal)")) %>%
    style_gt() %>%
    tab_source_note(source_note = md("Mann\u2013Whitney U per context (all DEGs and each intervention). \u201CNot estimable\u201D = no outlier STRs within that DEG context.")) %>%
    tab_source_note(source_note = legend_notes[[1]])
}

# ==========================================
# Panel D - Relative burden per genomic region
# ==========================================
gtD <- NULL
dr <- read_csv_safe(file.path(results_dir, "burden_region_mw.csv"))
if (!is.null(dr)) {
  pd <- data.table(
    Region = dr$region,
    `N (survivors/fatal)` = paste0(dr$n_survivors, "/", dr$n_fatal),
    `Survivors, median [IQR]` = dr$survivors_median_IQR,
    `Fatal, median [IQR]` = dr$fatal_median_IQR,
    `P (nominal)` = ifelse(dr$not_estimable, "Not estimable", fmt_pval(dr$p))
  )
  gtD <- pd %>%
    gt() %>%
    tab_header(title = md("**Supplementary Table S11 \u2014 Panel D. Relative burden per genomic region**")) %>%
    cols_label(Region = "Region", `N (survivors/fatal)` = "N (survivors/fatal)",
               `Survivors, median [IQR]` = "Survivors, median [IQR]",
               `Fatal, median [IQR]` = "Fatal, median [IQR]",
               `P (nominal)` = md("*P* (nominal)")) %>%
    style_gt() %>%
    tab_source_note(source_note = md("Mann\u2013Whitney U per genomic region. \u201CNot estimable\u201D = no outlier STRs detected in that region.")) %>%
    tab_source_note(source_note = legend_notes[[1]])
}

# ==========================================
# Save
# ==========================================
panels <- list(A = gtA, B = gtB, C = gtC, D = gtD)

for (nm in names(panels)) {
  if (is.null(panels[[nm]])) next
  f <- file.path(out_dir, sprintf("Supplementary_Table_S11_%s.html", nm))
  gtsave(panels[[nm]], f)
  cat(sprintf("  Saved: %s\n", f))
}

# Combined document
raw <- vapply(panels[!vapply(panels, is.null, logical(1))],
              function(g) as_raw_html(g), character(1))
if (length(raw) > 0) {
  combined <- paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>Supplementary Table S11</title></head><body>",
    paste(raw, collapse = "\n<hr>\n"),
    "</body></html>"
  )
  f_all <- file.path(out_dir, "Supplementary_Table_S11.html")
  writeLines(combined, f_all)
  cat(sprintf("  Saved: %s\n", f_all))
} else {
  cat("  [WARN] No panels generated (missing inputs).\n")
}

cat("\nDone.\n")
