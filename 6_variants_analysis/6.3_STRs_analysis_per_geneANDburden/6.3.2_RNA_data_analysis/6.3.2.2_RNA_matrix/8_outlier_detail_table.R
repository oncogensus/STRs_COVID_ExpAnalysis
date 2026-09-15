#!/usr/bin/env Rscript
# 8_outlier_detail_table.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Gera tabela publication-ready (gt HTML) com detalhes dos outliers
#   DBSCAN identificados no raincloud. Inclui localizacao genomica,
#   tamanhos de alelos, metricas DBSCAN e informacoes DEG.
#
# ENTRADAS (por argumentos de linha de comando)
#   --intervention-outliers  intervention_outliers.tsv (saida do step 1)
#   --out-dir                Diretorio de saida
#
# SAIDAS
#   outlier_detail_table.html
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(gt)
  library(scales)
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

path_outliers <- parse_arg("--intervention-outliers")
out_dir       <- parse_arg("--out-dir", ".")

if (is.null(path_outliers)) {
  stop("Argumento ausente: --intervention-outliers")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("--- Outlier Detail Table (Publication-Ready) ---\n")

dt <- fread(path_outliers, header = TRUE, sep = "\t")
cat(sprintf("  intervention_outliers.tsv: %d linhas\n", nrow(dt)))

# ==========================================
# 2. Rename DBSCAN columns (remove suffix)
# ==========================================
db_cols <- grep("_dbscan_global$", names(dt), value = TRUE)
if (length(db_cols) > 0) {
  short_names <- sub("_dbscan_global$", "", db_cols)
  setnames(dt, db_cols, short_names)
}

# ==========================================
# 3. Apply publication-quality intervention labels
# ==========================================
intervention_labels <- c(
  "GSE183533:COVID_vs_CONTROL"    = "Fatal COVID-19 vs. Controls",
  "GSE188847:COVID_vs_CONTROL"    = "Non-Survivors vs. Controls",
  "GSE157103:COVID_ICU_vs_NonICU" = "ICU vs. Non-Critical",
  "GSE157103:HFD45_ajustado_ICU"  = "ICU-Adjusted HFD45",
  "GSE188847:ICUVENT_vs_CONTROL"  = "IMV vs. Controls"
)

key <- paste0(dt$gse, ":", dt$intervention)
dt[, Comparison := fifelse(key %in% names(intervention_labels),
                           intervention_labels[key], intervention)]

# ==========================================
# 4. Prepare table data (aggregated per STR_ID)
# ==========================================
cat("\nPreparando dados da tabela...\n")

# Format FDR for display
fmt_fdr <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
         ifelse(x < 0.05, sprintf("%.3f", x),
                sprintf("%.2f", x)))
}
fmt_noise <- function(x) sprintf("%.1f%%", x * 100)

# Format noise ratio as percentage
fmt_noise <- function(x) {
  sprintf("%.1f%%", x * 100)
}

# Truncate outlier samples for display
fmt_samples <- function(x, max_chars = 50) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "-"
  ifelse(nchar(x) > max_chars,
         paste0(substr(x, 1, max_chars), "..."),
         x)
}

# Aggregate by STR_ID: one row per unique STR per GSE × Comparison
tbl <- dt[, .(
  GSE = first(gse),
  Comparison = first(Comparison),
  Gene = first(gene_name),
  Region = first(region),
  Chrom = first(chrom),
  Start = first(start),
  Motif = first(repeat_unit),
  Allele_1_Med = round(median(allele1_est), 1),
  Allele_1_Min = min(allele1_est),
  Allele_1_Max = max(allele1_est),
  Allele_2_Med = round(median(allele2_est), 1),
  Allele_2_Min = min(allele2_est),
  Allele_2_Max = max(allele2_est),
  Depth_Med = round(median(depth), 1),
  Depth_Min = min(depth),
  Depth_Max = max(depth),
  Clusters = first(n_clusters),
  `Noise %` = fmt_noise(first(noise_ratio)),
  `N Outliers` = first(n_outliers),
  logFC = first(round(logFC, 2)),
  FDR = fmt_fdr(first(FDR))
), by = .(STRs_ID, gse, Comparison)]

setnames(tbl, "gse", "GSE")

# Create combined group column for gt (avoids duplicate column error)
tbl[, Group := paste0(GSE, " | ", Comparison)]

# Sort by GSE, Comparison, Gene
tbl <- tbl[order(GSE, Comparison, Gene)]

# Drop GSE, Comparison, STRs_ID (now encoded in Group)
tbl[, c("GSE", "Comparison", "STRs_ID") := NULL]

cat(sprintf("  STRs_ID unicos: %d (de %d observacoes originais)\n", nrow(tbl), nrow(dt)))
cat(sprintf("  Tabela final: %d linhas\n", nrow(tbl)))

# ==========================================
# 5. Build gt table
# ==========================================
cat("\nGerando tabela gt...\n")

n_comparisons <- uniqueN(tbl$Comparison)
n_genes <- uniqueN(tbl$Gene)
n_strs <- nrow(tbl)

out_gt <- tbl %>%
  gt(groupname_col = "Group") %>%
  tab_header(
    title = md("**DBSCAN Outlier Loci in RNA-Seq DEGs**"),
    subtitle = sprintf("%d outlier STR loci across %d comparisons and %d genes",
                       n_strs, n_comparisons, n_genes)
  ) %>%
  tab_spanner(
    label = "Genomic Location",
    columns = c(Gene, Region, Chrom, Start, Motif)
  ) %>%
  tab_spanner(
    label = "Allele 1 (median/min/max)",
    columns = c(Allele_1_Med, Allele_1_Min, Allele_1_Max)
  ) %>%
  tab_spanner(
    label = "Allele 2 (median/min/max)",
    columns = c(Allele_2_Med, Allele_2_Min, Allele_2_Max)
  ) %>%
  tab_spanner(
    label = "Depth (median/min/max)",
    columns = c(Depth_Med, Depth_Min, Depth_Max)
  ) %>%
  tab_spanner(
    label = "DBSCAN Metrics",
    columns = c(Clusters, `Noise %`, `N Outliers`)
  ) %>%
  tab_spanner(
    label = "Differential Expression",
    columns = c(logFC, FDR)
  ) %>%
  cols_label(
    Gene = "Gene",
    Region = "Region",
    Chrom = "Chr",
    Start = "Start",
    Motif = "Motif",
    Allele_1_Med = "Med",
    Allele_1_Min = "Min",
    Allele_1_Max = "Max",
    Allele_2_Med = "Med",
    Allele_2_Min = "Min",
    Allele_2_Max = "Max",
    Depth_Med = "Med",
    Depth_Min = "Min",
    Depth_Max = "Max",
    Clusters = "Clusters",
    `Noise %` = "Noise",
    `N Outliers` = "N",
    logFC = "logFC",
    FDR = "FDR"
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(
      cell_text(weight = "bold", size = px(10)),
      cell_borders(sides = "bottom", weight = px(1.5), color = "grey60")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(10)),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(size = px(9)),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(9)),
    locations = cells_body(columns = Gene)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "grey95"),
      cell_text(weight = "bold", size = px(11))
    ),
    locations = cells_row_groups()
  ) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(9),
    heading.align = "left",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.width = px(1.5),
    table_body.border.bottom.color = "black",
    table_body.hlines.color = "grey92",
    table_body.hlines.width = px(0.5),
    table.border.left.width = px(0),
    table.border.right.width = px(0),
    data_row.padding = px(3),
    row_group.padding = px(6)
  ) %>%
  tab_source_note(
    source_note = md("*Allele sizes in repeat units. Depth = sequencing coverage at locus. Values shown as median (min–max) per STR locus.*")
  ) %>%
  tab_source_note(
    source_note = md("*DBSCAN clusters = number of genotype clusters detected. Noise = fraction of unclustered observations.*")
  ) %>%
  tab_source_note(
    source_note = md("*FDR = Benjamini-Hochberg adjusted p-value from DEG analysis.*")
  ) %>%
  tab_source_note(
    source_note = md("*Source: intervention_outliers.tsv (cross_intervention_STRs.py)*")
  )

# ==========================================
# 6. Save
# ==========================================
out_html <- file.path(out_dir, "outlier_detail_table.html")
gtsave(out_gt, out_html)
cat(sprintf("\nTabela salva em: %s\n", out_html))

cat("\nConcluido.\n")
