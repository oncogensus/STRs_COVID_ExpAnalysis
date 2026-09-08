#!/usr/bin/env Rscript
# plot_rna_summary.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Gera visualizacoes para o cruzamento RNA-Seq x STRs:
#     1) Raincloud plot: distribuicao de allele2_est por estudo (GSE),
#        colorido por group (case/control), com outliers DBSCAN.
#     2) Tabela de publicacao: por GSE, contagens e proporcoes de
#        outliers DBSCAN e sobreposicao de alelos.
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog       STRs_analysis_dataset.tsv (per-sample x STR)
#   --rna-gene-strs     rna_gene_strs.tsv (STRs que sao DEGs, com datasets)
#   --rna-outliers      rna_outlier_genes.tsv (outliers DBSCAN global)
#   --summary           rna_summary_by_study.tsv (resumo por estudo x gene)
#   --out-dir           Diretorio de saida
#
# SAIDAS
#   rna_raincloud_by_study.png   Raincloud plot (apenas outliers DBSCAN)
#   rna_publication_table.tsv    Tabela por GSE
#   rna_publication_table.html   Tabela gt formatada (HTML)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrain)
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

path_str_catalog  <- parse_arg("--str-catalog")
path_rna_gene_strs <- parse_arg("--rna-gene-strs")
path_rna_outliers <- parse_arg("--rna-outliers")
path_summary      <- parse_arg("--summary")
out_dir           <- parse_arg("--out-dir", ".")

missing <- c(
  "--str-catalog" = is.null(path_str_catalog),
  "--rna-gene-strs" = is.null(path_rna_gene_strs),
  "--rna-outliers" = is.null(path_rna_outliers),
  "--summary" = is.null(path_summary)
)
if (any(missing)) {
  stop("Argumentos ausentes: ", paste(names(missing[missing]), collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("Carregando STR catalog...\n")
str_cat <- fread(path_str_catalog, header = TRUE, sep = "\t")
cat(sprintf("  %d linhas (sample x STR), %d colunas\n", nrow(str_cat), ncol(str_cat)))

cat("Carregando rna_gene_strs.tsv...\n")
rna_gene_strs <- fread(path_rna_gene_strs, header = TRUE, sep = "\t")
cat(sprintf("  %d pares gene x STR\n", nrow(rna_gene_strs)))

cat("Carregando rna_outlier_genes.tsv...\n")
rna_outliers <- fread(path_rna_outliers, header = TRUE, sep = "\t")
cat(sprintf("  %d linhas (gene x STR x GSE)\n", nrow(rna_outliers)))

cat("Carregando rna_summary_by_study.tsv...\n")
rna_summary <- fread(path_summary, header = TRUE, sep = "\t")
cat(sprintf("  %d linhas (gene x estudo)\n", nrow(rna_summary)))

# ==========================================
# 2. Prepare data
# ==========================================
cat("\nPreparando dados...\n")

de_strs <- unique(rna_gene_strs$strs_id)
cat(sprintf("  STRs unicos em DEGs: %d\n", length(de_strs)))

str_deg <- str_cat[STRs_ID %in% de_strs]
cat(sprintf("  Linhas no catalogo para DEGs: %d\n", nrow(str_deg)))

if (!"group" %in% colnames(str_deg)) {
  stop("Coluna 'group' nao encontrada no STR catalog. Verifique 6.1_merge_datasets.r")
}
cat(sprintf("  Grupos encontrados: %s\n",
            paste(unique(str_deg$group), collapse = ", ")))

gse_map <- rna_gene_strs[, .(strs_id, datasets)]
gse_map[, datasets := trimws(datasets)]
gse_map <- gse_map[, .(gse = unlist(tstrsplit(datasets, ";", fixed = TRUE))),
                   by = strs_id]
gse_map <- gse_map[gse != ""]

str_deg <- merge(str_deg, gse_map, by.x = "STRs_ID", by.y = "strs_id",
                 allow.cartesian = TRUE)
cat(sprintf("  Apos merge com GSE: %d linhas\n", nrow(str_deg)))

outlier_strs <- unique(rna_outliers$strs_id)
str_deg[, is_outlier := STRs_ID %in% outlier_strs]
cat(sprintf("  STRs com outlier DBSCAN global: %d\n",
            sum(unique(str_deg[, .(STRs_ID, is_outlier)])$is_outlier)))

str_deg <- merge(str_deg,
                 rna_summary[, .(gse, gene, overlap_maior_alealo_grupos)],
                 by.x = c("gse", "gene_name"), by.y = c("gse", "gene"),
                 all.x = TRUE)

# ==========================================
# 3. Raincloud plot (scientific theme)
# ==========================================
cat("\nGerando raincloud plot...\n")

group_colors <- c("case" = "#E41A1C", "control" = "#377EB8")

str_deg_plot <- str_deg[is_outlier == TRUE & !is.na(group) & group != ""]

gene_counts <- str_deg_plot[, .N, by = gene_name]
valid_genes <- gene_counts[N >= 10, gene_name]
str_deg_plot <- str_deg_plot[gene_name %in% valid_genes]

str_deg_plot[, group := factor(group, levels = c("case", "control"))]

cat(sprintf("  Variantes para raincloud (apenas outliers DBSCAN global, >=10 obs): %d linhas, %d genes\n",
            nrow(str_deg_plot), length(unique(str_deg_plot$gene_name))))

p_rain <- ggplot(str_deg_plot,
                 aes(x = allele2_est, y = gene_name, fill = group)) +
  geom_rain(
    alpha = 0.6,
    point.args = list(size = 1.2, alpha = 0.5),
    boxplot.args = list(outlier.shape = NA, width = 0.2)
  ) +
  facet_wrap(~ gse, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 1, 5, 10, 50, 100, 500, 1000),
    labels = c("0", "1", "5", "10", "50", "100", "500", "1k"),
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Allele 2 length distributions in DEG STRs",
    subtitle = "Per GSE study, colored by case/control status (DBSCAN global outliers only)",
    x = "Allele 2 length (repeat units, log1p)",
    y = "Gene",
    caption = "Only DBSCAN global outliers shown"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "grey20"),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold",
                              margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, color = "grey40",
                                 margin = margin(b = 10)),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

out_png <- file.path(out_dir, "rna_raincloud_by_study.png")
raincloud_h <- min(40, max(6, length(unique(str_deg_plot$gene_name)) * 0.5 + 2))
ggsave(
  filename = out_png,
  plot = p_rain,
  width = 12,
  height = raincloud_h,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
cat(sprintf("Raincloud salvo em: %s\n", out_png))

# ==========================================
# 4. Publication table (per GSE)
# ==========================================
cat("\nGerando tabela de publicacao por GSE...\n")

# --- 4.1 Outlier counts per GSE x group (sample-level observations) ---
outlier_by_gse <- str_deg[is_outlier == TRUE & !is.na(group), .(
  n_outliers = .N,
  n_out_case = sum(group == "case"),
  n_out_control = sum(group == "control")
), by = gse]

# --- 4.2 Outliers with no overlap (sample-level observations) ---
out_no_overlap <- str_deg[is_outlier == TRUE & overlap_maior_alealo_grupos == "nao" &
                            !is.na(group), .(
  n_out_no_overlap = .N,
  n_out_no_overlap_case = sum(group == "case"),
  n_out_no_overlap_control = sum(group == "control")
), by = gse]

# --- 4.3 STR density per gene per GSE ---
rna_gene_strs_gse <- rna_gene_strs[, .(strs_id, gene, datasets)]
rna_gene_strs_gse[, datasets := trimws(datasets)]
rna_gene_strs_gse <- rna_gene_strs_gse[, .(gse = unlist(tstrsplit(datasets, ";", fixed = TRUE))),
                                        by = .(strs_id, gene)]
rna_gene_strs_gse <- rna_gene_strs_gse[gse != ""]
strs_per_gene <- rna_gene_strs_gse[, .(n_strs = .N), by = .(gse, gene)]
str_density <- strs_per_gene[, .(
  str_density_per_gene = round(sum(n_strs) / uniqueN(gene), 2)
), by = gse]

# --- 4.4 Summary counts from rna_summary ---
summary_gse <- rna_summary[, .(
  n_genes = uniqueN(gene),
  n_strs_total = sum(n_strs_identified),
  n_overlap_sim = sum(overlap_maior_alealo_grupos == "sim", na.rm = TRUE),
  n_overlap_nao = sum(overlap_maior_alealo_grupos == "nao", na.rm = TRUE),
  n_sem_dados = sum(overlap_maior_alealo_grupos == "sem_dados", na.rm = TRUE)
), by = gse]

# --- 4.5 Case/control counts from str_deg ---
group_counts <- str_deg[!is.na(group), .(
  n_case = sum(group == "case"),
  n_control = sum(group == "control")
), by = gse]

# --- 4.6 Median allele2_est per group per GSE ---
allele_medians <- str_deg[!is.na(group), .(
  median_allele_case = round(median(allele2_est[group == "case"], na.rm = TRUE), 2),
  median_allele_control = round(median(allele2_est[group == "control"], na.rm = TRUE), 2)
), by = gse]

# --- 4.7 Merge all ---
pub_table <- merge(summary_gse, str_density, by = "gse", all.x = TRUE)
pub_table <- merge(pub_table, group_counts, by = "gse", all.x = TRUE)
pub_table <- merge(pub_table, allele_medians, by = "gse", all.x = TRUE)
pub_table <- merge(pub_table, outlier_by_gse, by = "gse", all.x = TRUE)
pub_table <- merge(pub_table, out_no_overlap, by = "gse", all.x = TRUE)

# Fill NA with 0
for (col in c("n_outliers", "n_out_case", "n_out_control",
              "n_out_no_overlap", "n_out_no_overlap_case", "n_out_no_overlap_control")) {
  pub_table[is.na(get(col)), (col) := 0]
}

# --- 4.8 Format columns: abs (pct%) ---
fmt_pct <- function(x, total) {
  fifelse(total > 0,
          sprintf("%d (%.1f%%)", x, x / total * 100),
          "0 (0.0%)")
}

pub_table[, n_case_str := fmt_pct(n_case, n_case + n_control)]
pub_table[, n_control_str := fmt_pct(n_control, n_case + n_control)]
pub_table[, n_outliers_str := fmt_pct(n_outliers, n_strs_total)]
pub_table[, n_out_case_str := fmt_pct(n_out_case, n_strs_total)]
pub_table[, n_out_control_str := fmt_pct(n_out_control, n_strs_total)]
pub_table[, n_overlap_str := fmt_pct(n_overlap_sim, n_strs_total)]
pub_table[, n_out_no_overlap_str := fmt_pct(n_out_no_overlap, n_strs_total)]

# Order by n_outliers descending
pub_table <- pub_table[order(-n_outliers)]

cat(sprintf("  Tabela final: %d linhas (estudos GSE)\n", nrow(pub_table)))
cat(sprintf("  Estudos: %s\n", paste(pub_table$gse, collapse = ", ")))

# Save TSV
out_tsv <- file.path(out_dir, "rna_publication_table.tsv")
pub_tsv <- pub_table[, .(gse, n_genes, n_strs_total, str_density_per_gene,
                          n_case_str, n_control_str,
                          median_allele_case, median_allele_control,
                          n_outliers_str, n_out_case_str, n_out_control_str,
                          n_overlap_str, n_out_no_overlap_str)]
cat("  Primeiras 3 linhas:\n")
print(head(as.data.frame(pub_tsv), 3))
fwrite(pub_tsv, out_tsv, sep = "\t")
cat(sprintf("Tabela salva em: %s (%d estudos)\n", out_tsv, nrow(pub_tsv)))

# Save formatted gt table as HTML
pub_gt <- pub_tsv %>%
  gt() %>%
  tab_header(
    title = md("**RNA-Seq x STRs: Summary per GSE Study**"),
    subtitle = "STRs in DEG genes, DBSCAN outliers, and allele overlap by group"
  ) %>%
  cols_label(
    gse = "GSE",
    n_genes = "DEG Genes",
    n_strs_total = "STRs Total",
    str_density_per_gene = "STR Density/gene",
    n_case_str = "Case",
    n_control_str = "Control",
    median_allele_case = "Median Allele2 (case)",
    median_allele_control = "Median Allele2 (control)",
    n_outliers_str = "Outliers DBSCAN",
    n_out_case_str = "Outliers (case)",
    n_out_control_str = "Outliers (control)",
    n_overlap_str = "Overlap alelo maior",
    n_out_no_overlap_str = "Outliers sem overlap"
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_source_note(
    source_note = "Values: absolute count (percentage%). Overlap = same allele range between case/control."
  ) %>%
  tab_source_note(
    source_note = "Outliers sem overlap = outlier observations where allele ranges do not overlap between groups."
  ) %>%
  tab_source_note(
    source_note = "Source: cross_DEGs_STRs.py output"
  ) %>%
  opt_stylize(style = 3)

out_gt_html <- file.path(out_dir, "rna_publication_table.html")
gtsave(pub_gt, out_gt_html)
cat(sprintf("Tabela gt salva em: %s\n", out_gt_html))

cat("\nConcluido.\n")
