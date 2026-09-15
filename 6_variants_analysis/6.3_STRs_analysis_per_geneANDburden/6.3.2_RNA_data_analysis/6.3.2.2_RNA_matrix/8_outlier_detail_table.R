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

# Count unique outlier samples per group per comparison
dt[, sample_id := sub(";$", "", outlier_samples)]
n_per_group <- unique(dt[, .(gse, Comparison, group, sample_id)])[, .N, by = .(gse, Comparison, group)]
n_wide <- dcast(n_per_group, gse + Comparison ~ group, value.var = "N", fill = 0)
if ("case" %in% names(n_wide) && "control" %in% names(n_wide)) {
  setnames(n_wide, c("case", "control"), c("n_cases", "n_controls"))
}

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

# Prettify region labels
region_labels <- c(
  "three_prime_utr"  = "3\u2032 UTR",
  "five_prime_utr"   = "5\u2032 UTR",
  "promoter"         = "Promoter",
  "intron"           = "Intron",
  "exon"             = "Exon",
  "non_coding_exons" = "Non-coding exon",
  "intergenic"       = "Intergenic",
  "others"           = "Other"
)
dt[, region := fifelse(region %in% names(region_labels),
                       region_labels[region], region)]

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
# Lookup sample counts directly (avoids merge duplicate column issue)
n_wide_key <- paste0(n_wide$gse, "||", n_wide$Comparison)
names(n_wide_key) <- NULL
n_cases_vec <- setNames(n_wide$n_cases, n_wide_key)
n_controls_vec <- setNames(n_wide$n_controls, n_wide_key)
tbl_key <- paste0(tbl$GSE, "||", tbl$Comparison)
tbl[, n_cases := n_cases_vec[tbl_key]]
tbl[, n_controls := n_controls_vec[tbl_key]]
tbl[, n_cases := fifelse(is.na(n_cases), 0L, as.integer(n_cases))]
tbl[, n_controls := fifelse(is.na(n_controls), 0L, as.integer(n_controls))]
tbl[, Group := paste0(GSE, " | ", Comparison, " | Cases: ", n_cases, " / Controls: ", n_controls)]

# Sort by GSE, Comparison, Gene
tbl <- tbl[order(GSE, Comparison, Gene)]

# Drop GSE, Comparison, STRs_ID (now encoded in Group)
tbl[, c("GSE", "Comparison", "STRs_ID", "n_cases", "n_controls") := NULL]

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
    label = md("Allele 1 (median/min/max)^a^"),
    columns = c(Allele_1_Med, Allele_1_Min, Allele_1_Max)
  ) %>%
  tab_spanner(
    label = md("Allele 2 (median/min/max)^a^"),
    columns = c(Allele_2_Med, Allele_2_Min, Allele_2_Max)
  ) %>%
  tab_spanner(
    label = md("Depth (median/min/max)^a^"),
    columns = c(Depth_Med, Depth_Min, Depth_Max)
  ) %>%
  tab_spanner(
    label = md("DBSCAN Metrics^b^"),
    columns = c(Clusters, `Noise %`, `N Outliers`)
  ) %>%
  tab_spanner(
    label = md("Differential Expression^c^"),
    columns = c(logFC, FDR)
  ) %>%
  cols_label(
    Gene = "Gene",
    Region = "Region",
    Chrom = "Chr",
    Start = "Start",
    Motif = "Motif",
    Allele_1_Med = md("Med^a^"),
    Allele_1_Min = md("Min^a^"),
    Allele_1_Max = md("Max^a^"),
    Allele_2_Med = md("Med^a^"),
    Allele_2_Min = md("Min^a^"),
    Allele_2_Max = md("Max^a^"),
    Depth_Med = md("Med^a^"),
    Depth_Min = md("Min^a^"),
    Depth_Max = md("Max^a^"),
    Clusters = md("Clusters^b^"),
    `Noise %` = md("Noise^b^"),
    `N Outliers` = md("N^b^"),
    logFC = md("logFC^c^"),
    FDR = md("FDR^c^")
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
    source_note = md("*(a)* Allele sizes in repeat units. Depth = sequencing coverage at locus. Values represent median (min\u2013max) per locus across outlier samples.")
  ) %>%
  tab_source_note(
    source_note = md("*(b)* DBSCAN clustering metrics. Noise = fraction of unclustered observations (cluster 0). Outlier loci = STRs with \u22651 genotypic cluster detected.")
  ) %>%
  tab_source_note(
    source_note = md("*(c)* Log~2~ fold-change and Benjamini\u2013Hochberg adjusted *P*-value from differential expression analysis (DESeq2).")
  )

# ==========================================
# 6. Save
# ==========================================
out_html <- file.path(out_dir, "outlier_detail_table.html")
gtsave(out_gt, out_html)
cat(sprintf("\nTabela salva em: %s\n", out_html))

cat("\nConcluido.\n")
