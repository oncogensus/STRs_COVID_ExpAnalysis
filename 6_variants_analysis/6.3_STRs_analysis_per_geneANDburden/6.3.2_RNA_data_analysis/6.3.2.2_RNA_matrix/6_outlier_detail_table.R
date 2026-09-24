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
path_others_csv <- parse_arg("--others-csv")
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
  "GSE188847:COVID_vs_CONTROL"    = "Fatal COVID-19 vs. Controls",
  "GSE157103:COVID_ICU_vs_NonICU" = "ICU vs. Non-Critical",
  "GSE157103:HFD45_ajustado_ICU"  = "ICU-Adjusted HFD45",
  "GSE188847:ICUVENT_vs_CONTROL"  = "IMV vs. Controls"
)

key <- paste0(dt$gse, ":", dt$intervention)
dt[, Comparison := fifelse(key %in% names(intervention_labels),
                           intervention_labels[key], intervention)]

# Count outlier samples per group per variant (each row = 1 sample observation)
n_per_group <- dt[, .N, by = .(STRs_ID, gse, Comparison, group)]
n_wide <- dcast(n_per_group, STRs_ID + gse + Comparison ~ group, value.var = "N", fill = 0)
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
  "others"           = "Non-Coding elements"
)
dt[, region := fifelse(region %in% names(region_labels),
                       region_labels[region], region)]

# Lookup biotype for 'others' regions from annotation TSV
if (!is.null(path_others_csv) && file.exists(path_others_csv)) {
  others_dt <- fread(path_others_csv, select = c("STRs_ID", "gene_biotype"))
  others_dt <- others_dt[gene_biotype != "" & !is.na(gene_biotype)]
  others_dt <- unique(others_dt, by = "STRs_ID")
  biotype_map <- setNames(others_dt$gene_biotype, others_dt$STRs_ID)
  other_rows <- which(dt$region == "Other")
  if (length(other_rows) > 0) {
    matched <- biotype_map[dt$STRs_ID[other_rows]]
    n_reclass <- sum(!is.na(matched))
    dt$region[other_rows] <- ifelse(is.na(matched), "Non-coding", matched)
    cat(sprintf("  Biotype lookup: %d/%d variantes 'others' reclassificados\n",
                n_reclass, length(other_rows)))
  }
}

# Aggregate by STR_ID: one row per unique STR per GSE × Comparison
tbl <- dt[, .(
  Gene = first(gene_name),
  Region = first(region),
  Chrom = first(chrom),
  Start = first(start),
  Motif = first(repeat_unit),
  Depth_Med = round(median(depth), 1),
  Depth_Min = min(depth),
  Depth_Max = max(depth),
  Clusters = first(n_clusters),
  `Noise %` = fmt_noise(first(noise_ratio)),
  `N Outliers` = first(n_outliers),
  logFC = first(round(logFC, 2)),
  FDR = fmt_fdr(first(FDR))
), by = .(STRs_ID, gse, Comparison)]

# Allele 2 sizes separately per group (case / control)
aux_case <- dt[group == "case", .(
  Allele_2_Case_Med = round(median(allele2_est), 1),
  Allele_2_Case_Min = min(allele2_est),
  Allele_2_Case_Max = max(allele2_est)
), by = .(STRs_ID, gse, Comparison)]

aux_control <- dt[group == "control", .(
  Allele_2_Control_Med = round(median(allele2_est), 1),
  Allele_2_Control_Min = min(allele2_est),
  Allele_2_Control_Max = max(allele2_est)
), by = .(STRs_ID, gse, Comparison)]

tbl <- merge(tbl, aux_case, by = c("STRs_ID", "gse", "Comparison"), all.x = TRUE)
tbl <- merge(tbl, aux_control, by = c("STRs_ID", "gse", "Comparison"), all.x = TRUE)

setnames(tbl, "gse", "GSE")

# Create combined group column for gt (avoids duplicate column error)
# Lookup sample counts per STR_ID (avoids merge duplicate column issue)
n_wide_key <- paste0(n_wide$STRs_ID, "||", n_wide$gse, "||", n_wide$Comparison)
n_cases_vec <- setNames(n_wide$n_cases, n_wide_key)
n_controls_vec <- setNames(n_wide$n_controls, n_wide_key)
tbl_key <- paste0(tbl$STRs_ID, "||", tbl$GSE, "||", tbl$Comparison)
tbl[, n_cases := as.integer(n_cases_vec[tbl_key])]
tbl[, n_controls := as.integer(n_controls_vec[tbl_key])]
tbl[, n_cases := fifelse(is.na(n_cases), 0L, n_cases)]
tbl[, n_controls := fifelse(is.na(n_controls), 0L, n_controls)]
tbl[, Group := paste0(GSE, " | ", Comparison)]

# Rename count columns
setnames(tbl, c("n_cases", "n_controls"), c("Cases", "Controls"))

# Sort by GSE, Comparison, Gene
tbl <- tbl[order(GSE, Comparison, Gene)]

# Compute stats before dropping columns
n_comparisons <- uniqueN(tbl$Comparison)
n_genes <- uniqueN(tbl$Gene)
n_strs <- nrow(tbl)

# Drop auxiliary columns — keep only columns that go into the table
keep_cols <- c("Group", "Gene", "Region", "Chrom", "Start", "Motif",
               paste0("Allele_2_Case_", c("Med", "Min", "Max")),
               paste0("Allele_2_Control_", c("Med", "Min", "Max")),
               paste0("Depth_", c("Med", "Min", "Max")),
               "Clusters", "Noise %", "N Outliers", "logFC", "FDR",
               "Cases", "Controls")
tbl <- tbl[, ..keep_cols]

cat(sprintf("  STRs_ID unicos: %d (de %d observacoes originais)\n", n_strs, nrow(dt)))
cat(sprintf("  Tabela final: %d linhas\n", nrow(tbl)))

# ==========================================
# 5. Build gt table
# ==========================================
cat("\nGerando tabela gt...\n")

out_gt <- tbl %>%
  gt(groupname_col = "Group") %>%
  tab_header(
    title = md("**DBSCAN Outliers Loci in RNA-Seq DEGs**"),
  ) %>%
  tab_spanner(
    label = "Genomic Location",
    columns = c(Gene, Region, Chrom, Start, Motif)
  ) %>%
  tab_spanner(
    label = md("Allele 2 - Case"),
    columns = c(Allele_2_Case_Med, Allele_2_Case_Min, Allele_2_Case_Max)
  ) %>%
  tab_spanner(
    label = md("Allele 2 - Control"),
    columns = c(Allele_2_Control_Med, Allele_2_Control_Min, Allele_2_Control_Max)
  ) %>%
  tab_spanner(
    label = md("Depth^a^"),
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
  tab_spanner(
    label = "Sample Counts",
    columns = c(Cases, Controls)
  ) %>%
  cols_label(
    Gene = "Gene",
    Region = "Region",
    Chrom = "Chr",
    Start = "Start",
    Motif = "Motif",
    Allele_2_Case_Med = md("Med^a^"),
    Allele_2_Case_Min = md("Min^a^"),
    Allele_2_Case_Max = md("Max^a^"),
    Allele_2_Control_Med = md("Med^a^"),
    Allele_2_Control_Min = md("Min^a^"),
    Allele_2_Control_Max = md("Max^a^"),
    Depth_Med = md("Med^a^"),
    Depth_Min = md("Min^a^"),
    Depth_Max = md("Max^a^"),
    Clusters = md("Clusters^b^"),
    `Noise %` = md("Noise^b^"),
    `N Outliers` = md("N^b^"),
    logFC = md("logFC^c^"),
    FDR = md("FDR^c^"),
    Cases = "Cases",
    Controls = "Controls"
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "grey95"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
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
  ) %>%
  tab_source_note(
    source_note = md("*(a)* Allele 2 sizes in repeat units per group (case/control). Depth = sequencing coverage at locus. Values represent median (min\u2013max) per locus across outlier samples.")
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
