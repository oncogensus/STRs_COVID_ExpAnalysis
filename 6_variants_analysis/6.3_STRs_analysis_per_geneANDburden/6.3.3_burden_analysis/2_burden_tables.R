#!/usr/bin/env Rscript
# 2_burden_tables.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Builds publication-ready supplementary tables (gt HTML) from the burden
#   analysis outputs (step 1_burden_analysis.R). Four panels:
#     A - Global relative burden (Mann-Whitney)
#     B - Firth logistic regression of global relative burden
#     C - Relative burden within DEGs (global + per intervention)
#     D - Relative burden per genomic region
#   Panels report nominal P; a Benjamini-Hochberg FDR is shown for B, C and D
#   (Panel A is a single comparison, nominal only). Rows without variance in
#   both groups (IQR span = 0) are omitted. Style matches the other
#   supplementary tables.
#
# INPUTS (via command-line arguments)
#   --results-dir   Directory with burden_*.csv (default: results)
#   --out-dir       Output directory for tables (default: <results-dir>/tables)
#
# OUTPUTS
#   burden_analysis.html   All four panels in one document
#   burden_global.html     Panel A
#   burden_firth.html      Panel B
#   burden_deg.html        Panel C
#   burden_region.html     Panel D
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
# Formatters / helpers
# ==========================================
fmt_pval <- function(x) {
  ifelse(is.na(x), "-",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                ifelse(x < 0.05, sprintf("%.3f", x), sprintf("%.2f", x))))
}

fmt_num <- function(x, d = 2) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))

fmt_or_ci <- function(or, lo, hi, d = 2) {
  ifelse(is.na(or), "-", sprintf("%s (%s\u2013%s)", fmt_num(or, d), fmt_num(lo, d), fmt_num(hi, d)))
}

# Significance: ** FDR < 0.05; * nominal p < 0.05 (uncorrected).
sig_stars <- function(p, p_adj) {
  ifelse(!is.na(p_adj) & p_adj < 0.05, "**",
         ifelse(!is.na(p) & p < 0.05, "*", ""))
}

# Split "median [IQR]" strings into two vectors.
split_med <- function(s) trimws(sub("\\s*\\[.*$", "", s))
split_iqr <- function(s) {
  out <- trimws(sub("^[^[]*", "", s))
  ifelse(is.na(s) | !grepl("\\[", s), "", out)
}

# Width of the IQR interval (upper - lower). 0 = no variance (proxy), NA = missing.
# Used to keep only rows with variance in both groups (script-2-only heuristic).
iqr_span <- function(iqr) {
  inner <- gsub("\\[|\\]", "", iqr)
  vapply(strsplit(inner, "-", fixed = TRUE), function(v) {
    v <- suppressWarnings(as.numeric(v))
    v <- v[!is.na(v)]
    if (length(v) >= 2) abs(v[length(v)] - v[1]) else NA_real_
  }, numeric(1))
}

read_csv_safe <- function(f) {
  if (file.exists(f)) return(fread(f))
  cat(sprintf("  [WARN] missing: %s\n", f))
  NULL
}

# Common publication style (matches the other supplementary tables).
style_gt <- function(g) {
  g %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_spanners()
    ) %>%
    tab_options(
      table.font.size = px(13),
      table.width = pct(100),
      heading.align = "center",
      heading.title.font.size = px(16),
      heading.subtitle.font.size = px(12),
      column_labels.font.weight = "bold",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      table_body.border.bottom.style = "solid",
      row_group.padding = px(6),
      data_row.padding = px(5)
    )
}

bold_sig <- function(g) {
  tab_style(
    g,
    style = cell_text(weight = "bold", size = px(20)),
    locations = cells_body(columns = Sig, rows = Sig != "")
  )
}

note_sig <- md("Significance: <sup>**</sup> *P* < 0.05 (FDR); <sup>*</sup> *P* < 0.05 (nominal, uncorrected).")
note_legend <- md("*DEG*, differentially expressed gene; *IQR*, interquartile range; *BH*, Benjamini\u2013Hochberg; *OR*, odds ratio; *CI*, confidence interval; *PC1*, first ancestry principal component.")

cat("--- Supplementary burden tables ---\n")

# ==========================================
# Panel A - Global relative burden (Mann-Whitney)
# ==========================================
gtA <- NULL
gmw <- read_csv_safe(file.path(results_dir, "burden_global_mw.csv"))
if (!is.null(gmw)) {
  med <- c(gmw$survivors_median_IQR[1], gmw$fatal_median_IQR[1])
  spans <- iqr_span(split_iqr(med))
  if (all(!is.na(spans) & spans > 0)) {
    pa <- data.table(
      Group = c("Survivors", "Fatal COVID-19 cases"),
      N = c(gmw$n_survivors[1], gmw$n_fatal[1]),
      Median = split_med(med),
      IQR = split_iqr(med),
      Sig = c(sig_stars(gmw$p[1], NA_real_), ""),
      `P (nominal)` = c(fmt_pval(gmw$p[1]), "")
    )
    gtA <- pa %>%
      gt() %>%
      tab_header(title = md("**Global relative burden of outlier STRs between fatal cases and survivors**")) %>%
      tab_spanner(label = "Relative burden", columns = c(Median, IQR)) %>%
      cols_label(Group = "Group", N = "N", Median = "Median", IQR = "IQR", Sig = "",
                 `P (nominal)` = md("*P* (nominal)")) %>%
      style_gt() %>%
      bold_sig() %>%
      tab_source_note(source_note = md("Mann\u2013Whitney U test (single comparison; no multiple-testing correction).")) %>%
      tab_source_note(source_note = note_legend)
  } else {
    cat("  [WARN] Panel A skipped: no variance in one group.\n")
  }
}

# ==========================================
# Panel B - Firth logistic regression
# ==========================================
gtB <- NULL
fb <- read_csv_safe(file.path(results_dir, "burden_firth.csv"))
if (!is.null(fb)) {
  p_adj <- if ("p_adj" %in% names(fb)) fb$p_adj else rep(NA_real_, nrow(fb))
  pb <- data.table(
    Predictor = fb$predictor,
    `OR (95% CI)` = fmt_or_ci(fb$OR, fb$CI_low, fb$CI_high),
    Sig = sig_stars(fb$p, p_adj),
    `P (nominal)` = fmt_pval(fb$p),
    FDR = fmt_pval(p_adj)
  )
  gtB <- pb %>%
    gt() %>%
    tab_header(title = md("**Firth logistic regression of global relative burden on COVID-19 fatality**")) %>%
    cols_label(Predictor = "Predictor", `OR (95% CI)` = md("OR (95% CI)"), Sig = "",
               `P (nominal)` = md("*P* (nominal)"), FDR = "FDR") %>%
    style_gt() %>%
    bold_sig() %>%
    tab_source_note(source_note = md("Outcome: fatal COVID-19 (1) vs survivor (0). Relative burden per 1 percentage point; age per year; sex (male vs female); PC1 per standard deviation. Firth's penalized likelihood. FDR: BH across the four non-intercept predictors.")) %>%
    tab_source_note(source_note = note_legend)
}

# ==========================================
# Panel C - Relative burden within DEGs
# ==========================================
gtC <- NULL
dc <- read_csv_safe(file.path(results_dir, "burden_deg_mw.csv"))
if (!is.null(dc)) {
  p_adj <- if ("p_adj" %in% names(dc)) dc$p_adj else p.adjust(dc$p, method = "BH", n = sum(!is.na(dc$p)))
  pc <- data.table(
    Context = dc$context,
    `N (survivors/fatal)` = paste0(dc$n_survivors, "/", dc$n_fatal),
    s_med = split_med(dc$survivors_median_IQR),
    s_iqr = split_iqr(dc$survivors_median_IQR),
    f_med = split_med(dc$fatal_median_IQR),
    f_iqr = split_iqr(dc$fatal_median_IQR),
    Sig = sig_stars(dc$p, p_adj),
    `P (nominal)` = ifelse(is.na(dc$p), "Not estimable", fmt_pval(dc$p)),
    FDR = ifelse(is.na(dc$p), "-", fmt_pval(p_adj))
  )
  # Keep only contexts with variance in both groups (IQR span > 0).
  keep <- iqr_span(pc$s_iqr) > 0 & iqr_span(pc$f_iqr) > 0
  keep[is.na(keep)] <- FALSE
  pc <- pc[keep]
  gtC <- pc %>%
    gt() %>%
    tab_header(title = md("**Relative burden within DEGs by intervention**")) %>%
    tab_spanner(label = "Survivors", columns = c(s_med, s_iqr)) %>%
    tab_spanner(label = "Fatal COVID-19 cases", columns = c(f_med, f_iqr)) %>%
    cols_label(Context = "Context", `N (survivors/fatal)` = "N (survivors/fatal)",
               s_med = "Median", s_iqr = "IQR", f_med = "Median", f_iqr = "IQR",
               Sig = "", `P (nominal)` = md("*P* (nominal)"), FDR = "FDR") %>%
    style_gt() %>%
    bold_sig() %>%
    tab_source_note(source_note = md("Mann\u2013Whitney U per context (all DEGs and each intervention). \u201CNot estimable\u201D = no outlier STRs within that DEG context. Contexts without variance in both groups (IQR span = 0) were omitted. FDR: BH across all DEG contexts.")) %>%
    tab_source_note(source_note = note_sig)
}

# ==========================================
# Panel D - Relative burden per genomic region
# ==========================================
gtD <- NULL
dr <- read_csv_safe(file.path(results_dir, "burden_region_mw.csv"))
if (!is.null(dr)) {
  p_adj <- if ("p_adj" %in% names(dr)) dr$p_adj else p.adjust(dr$p, method = "BH", n = sum(!is.na(dr$p)))
  pd <- data.table(
    Region = dr$region,
    `N (survivors/fatal)` = paste0(dr$n_survivors, "/", dr$n_fatal),
    s_med = split_med(dr$survivors_median_IQR),
    s_iqr = split_iqr(dr$survivors_median_IQR),
    f_med = split_med(dr$fatal_median_IQR),
    f_iqr = split_iqr(dr$fatal_median_IQR),
    Sig = ifelse(dr$not_estimable, "", sig_stars(dr$p, p_adj)),
    `P (nominal)` = ifelse(dr$not_estimable, "Not estimable", fmt_pval(dr$p)),
    FDR = ifelse(dr$not_estimable, "-", fmt_pval(p_adj))
  )
  # Keep only regions with variance in both groups (IQR span > 0).
  keep <- iqr_span(pd$s_iqr) > 0 & iqr_span(pd$f_iqr) > 0
  keep[is.na(keep)] <- FALSE
  pd <- pd[keep]
  gtD <- pd %>%
    gt() %>%
    tab_header(title = md("**Relative burden by genomic region**")) %>%
    tab_spanner(label = "Survivors", columns = c(s_med, s_iqr)) %>%
    tab_spanner(label = "Fatal COVID-19 cases", columns = c(f_med, f_iqr)) %>%
    cols_label(Region = "Region", `N (survivors/fatal)` = "N (survivors/fatal)",
               s_med = "Median", s_iqr = "IQR", f_med = "Median", f_iqr = "IQR",
               Sig = "", `P (nominal)` = md("*P* (nominal)"), FDR = "FDR") %>%
    style_gt() %>%
    bold_sig() %>%
    tab_source_note(source_note = md("Mann\u2013Whitney U per genomic region. \u201CNot estimable\u201D = no outlier STRs detected in that region. Regions without variance in both groups (IQR span = 0) were omitted. FDR: BH across estimable regions.")) %>%
    tab_source_note(source_note = note_sig)
}

# ==========================================
# Save
# ==========================================
panels <- list(
  A = list(gt = gtA, file = "burden_global.html"),
  B = list(gt = gtB, file = "burden_firth.html"),
  C = list(gt = gtC, file = "burden_deg.html"),
  D = list(gt = gtD, file = "burden_region.html")
)

for (nm in names(panels)) {
  if (is.null(panels[[nm]]$gt)) next
  f <- file.path(out_dir, panels[[nm]]$file)
  gtsave(panels[[nm]]$gt, f)
  cat(sprintf("  Saved: %s\n", f))
}

# Combined document
raw <- vapply(panels, function(p) if (is.null(p$gt)) NA_character_ else as_raw_html(p$gt), character(1))
raw <- raw[!is.na(raw)]
if (length(raw) > 0) {
  combined <- paste0(
    "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
    "<title>Burden analyses of outlier STRs and COVID-19 fatality</title>",
    "<style>body{background-color:white;}</style></head><body>",
    paste(raw, collapse = "\n<hr>\n"),
    "</body></html>"
  )
  f_all <- file.path(out_dir, "burden_analysis.html")
  writeLines(combined, f_all)
  cat(sprintf("  Saved: %s\n", f_all))
} else {
  cat("  [WARN] No panels generated (missing inputs).\n")
}

cat("\nDone.\n")
