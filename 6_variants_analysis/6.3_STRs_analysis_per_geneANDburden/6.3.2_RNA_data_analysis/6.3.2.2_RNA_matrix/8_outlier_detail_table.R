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
# 4. Prepare table data
# ==========================================
cat("\nPreparando dados da tabela...\n")

# Format FDR for display
fmt_fdr <- function(x) {
  fifelse(x < 0.001,
          formatC(x, format = "e", digits = 1),
          sprintf("%.3f", x))
}

# Format noise ratio as percentage
fmt_noise <- function(x) {
  sprintf("%.1f%%", x * 100)
}

# Truncate outlier samples for display
fmt_samples <- function(x, max_chars = 40) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "-"
  ifelse(nchar(x) > max_chars,
         paste0(substr(x, 1, max_chars), "..."),
         x)
}

# Build table
tbl <- dt[, .(
  Gene = gene_name,
  STRs_ID = STRs_ID,
  Region = region,
  Chrom = chrom,
  Start = format(start, big.mark = ","),
  End = format(end, big.mark = ","),
  Motif = repeat_unit,
  Allele_1 = allele1_est,
  Allele_2 = allele2_est,
  Depth = depth,
  Clusters = n_clusters,
  `Noise %` = fmt_noise(noise_ratio),
  `N Outliers` = n_outliers,
  `Outlier Samples` = fmt_samples(outlier_samples),
  logFC = round(logFC, 2),
  FDR = fmt_fdr(FDR),
  Direction = Direction,
  GSE = gse,
  Comparison = Comparison
)]

# Sort by FDR ascending (most significant first)
tbl <- tbl[order(Gene, Comparison)]

cat(sprintf("  Tabela final: %d linhas\n", nrow(tbl)))

# ==========================================
# 5. Build gt table
# ==========================================
cat("\nGerando tabela gt...\n")

n_comparisons <- uniqueN(tbl$Comparison)
n_genes <- uniqueN(tbl$Gene)

out_gt <- tbl %>%
  gt(groupname_col = "Comparison") %>%
  tab_header(
    title = md("**DBSCAN Outlier Loci in RNA-Seq DEGs**"),
    subtitle = sprintf("%d outlier loci across %d comparisons and %d genes",
                       nrow(tbl), n_comparisons, n_genes)
  ) %>%
  tab_spanner(
    label = "Genomic Location",
    columns = c(Gene, STRs_ID, Region, Chrom, Start, End, Motif)
  ) %>%
  tab_spanner(
    label = "Allele Estimates",
    columns = c(Allele_1, Allele_2, Depth)
  ) %>%
  tab_spanner(
    label = "DBSCAN Metrics",
    columns = c(Clusters, `Noise %`, `N Outliers`, `Outlier Samples`)
  ) %>%
  tab_spanner(
    label = "Differential Expression",
    columns = c(logFC, FDR, Direction)
  ) %>%
  cols_label(
    Gene = "Gene",
    STRs_ID = "STRs ID",
    Region = "Region",
    Chrom = "Chr",
    Start = "Start",
    End = "End",
    Motif = "Motif",
    Allele_1 = "Allele 1",
    Allele_2 = "Allele 2",
    Depth = "Depth",
    Clusters = "Clusters",
    `Noise %` = "Noise",
    `N Outliers` = "N",
    `Outlier Samples` = "Samples",
    logFC = "logFC",
    FDR = "FDR",
    Direction = "Dir."
  ) %>%
  cols_hide(columns = c(GSE)) %>%
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
    style = cell_text(style = "italic", size = px(9)),
    locations = cells_body(columns = STRs_ID)
  ) %>%
  opt_row_grouping(
    columns = "Comparison",
    row_group_style = list(
      cell_fill(color = "grey95"),
      cell_text(weight = "bold", size = px(11))
    )
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
    row_group.padding = px(6),
    source_note.font.size = px(8)
  ) %>%
  tab_source_note(
    source_note = md("*Allele sizes in repeat units. Depth = sequencing coverage at locus.*")
  ) %>%
  tab_source_note(
    source_note = md("*DBSCAN clusters = number of genotype clusters detected. Noise = fraction of unclustered observations.*")
  ) %>%
  tab_source_note(
    source_note = md("*FDR = Benjamini-Hochberg adjusted p-value from DEG analysis. Dir. = expression direction in severe COVID-19.*")
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
