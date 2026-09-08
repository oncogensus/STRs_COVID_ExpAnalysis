#!/usr/bin/env Rscript
# plot_rna_summary.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Gera visualizacoes para o cruzamento RNA-Seq x STRs:
#     1) Ridgeline plot: distribuicao de allele2_est por estudo (GSE),
#        colorido por group (case/control), com outliers DBSCAN marcados.
#     2) Tabela de publicacao: por gene, contagens e proporcoes de
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
#   rna_ridgeline_by_study.png   Ridgeline plot
#   rna_publication_table.tsv    Tabela por gene
#   rna_publication_table.png    Tabela gt formatada
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
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
# 2. Prepare data for ridgeline
# ==========================================
cat("\nPreparando dados para ridgeline...\n")

# Get unique STRs that are DEGs
de_strs <- unique(rna_gene_strs$strs_id)
cat(sprintf("  STRs unicos em DEGs: %d\n", length(de_strs)))

# Filter str_cat to only DEG STRs
str_deg <- str_cat[STRs_ID %in% de_strs]
cat(sprintf("  Linhas no catalogo para DEGs: %d\n", nrow(str_deg)))

# Check group column
if (!"group" %in% colnames(str_deg)) {
  stop("Coluna 'group' nao encontrada no STR catalog. Verifique 6.1_merge_datasets.r")
}
cat(sprintf("  Grupos encontrados: %s\n",
            paste(unique(str_deg$group), collapse = ", ")))

# Create GSE mapping from rna_gene_strs (datasets column has semicolon-separated GSEs)
gse_map <- rna_gene_strs[, .(strs_id, datasets)]
gse_map[, datasets := trimws(datasets)]
gse_map <- gse_map[, .(gse = unlist(tstrsplit(datasets, ";", fixed = TRUE))),
                   by = strs_id]
gse_map <- gse_map[gse != ""]

# Add GSE info to str_deg
str_deg <- merge(str_deg, gse_map, by.x = "STRs_ID", by.y = "strs_id",
                 allow.cartesian = TRUE)
cat(sprintf("  Apos merge com GSE: %d linhas\n", nrow(str_deg)))

# Add outlier status from rna_outliers
outlier_strs <- unique(rna_outliers$strs_id)
str_deg[, is_outlier := STRs_ID %in% outlier_strs]
cat(sprintf("  STRs com outlier DBSCAN global: %d\n", sum(unique(str_deg[, .(STRs_ID, is_outlier)])$is_outlier)))

# ==========================================
# 3. Ridgeline plot
# ==========================================
cat("\nGerando ridgeline plot...\n")

# Color palette: case = red, control = blue
group_colors <- c(
  "case"   = "#E41A1C",
  "control" = "#377EB8"
)

# Filter: only groups with data
str_deg_plot <- str_deg[!is.na(group) & group != ""]
str_deg_plot[, group := factor(group, levels = c("case", "control"))]

# Facet by GSE, y-axis = gene_name
p_ridge <- ggplot(str_deg_plot,
                  aes(x = allele2_est, y = gene_name, fill = group)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    rel_min_height = 0.005,
    color = "white",
    bandwidth = 1.5
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
    subtitle = "Per GSE study, colored by case/control status",
    x = "Allele 2 length (repeat units, log1p)",
    y = "Gene"
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(hjust = 1, size = 12),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(0.8, "lines")
  )

# Add outlier points if any
if (any(str_deg_plot$is_outlier)) {
  outliers_only <- str_deg_plot[is_outlier == TRUE]
  p_ridge <- p_ridge +
    geom_point(
      data = outliers_only,
      aes(x = allele2_est, y = gene_name),
      inherit.aes = FALSE,
      shape = 24, fill = "black", color = "white",
      size = 1.5, stroke = 0.5,
      position = position_nudge(y = 0.15)
    ) +
    labs(caption = "Black triangles = DBSCAN global outliers")
}

out_png <- file.path(out_dir, "rna_ridgeline_by_study.png")
ridgeline_h <- min(40, max(6, length(de_strs) * 0.4 + 2))
ggsave(
  filename = out_png,
  plot = p_ridge,
  width = 12,
  height = ridgeline_h,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
cat(sprintf("Ridgeline salvo em: %s\n", out_png))

# ==========================================
# 4. Publication table (per gene)
# ==========================================
cat("\nGerando tabela de publicacao...\n")

# For each gene, count STRs, outliers, overlap status
pub_table <- rna_summary[, .(
  n_strs_identified = sum(n_strs_identified),
  n_strs_outliers = sum(n_strs_identified_outliers),
  n_overlap_sim = sum(overlap_maior_alealo_grupos == "sim", na.rm = TRUE),
  n_overlap_nao = sum(overlap_maior_alealo_grupos == "nao", na.rm = TRUE),
  n_sem_dados = sum(overlap_maior_alealo_grupos == "sem_dados", na.rm = TRUE),
  gse = paste(unique(gse), collapse = ";")
), by = gene]

pub_table[, prop_outliers := fifelse(n_strs_identified > 0,
                                      sprintf("%.1f", n_strs_outliers / n_strs_identified * 100),
                                      "0.0")]
pub_table[, prop_overlap := fifelse(n_strs_identified > 0,
                                     sprintf("%.1f", n_overlap_sim / n_strs_identified * 100),
                                     "0.0")]
pub_table[, prop_sem_overlap := fifelse(n_strs_identified > 0,
                                         sprintf("%.1f", n_overlap_nao / n_strs_identified * 100),
                                         "0.0")]

# Order by n_strs_outliers descending
pub_table <- pub_table[order(-n_strs_outliers)]

# Save TSV
out_tsv <- file.path(out_dir, "rna_publication_table.tsv")
fwrite(pub_table, out_tsv, sep = "\t")
cat(sprintf("Tabela salva em: %s (%d genes)\n", out_tsv, nrow(pub_table)))

# Save formatted gt table as PNG
pub_gt <- pub_table %>%
  select(gene, gse, n_strs_identified, n_strs_outliers, prop_outliers,
         n_overlap_sim, n_sem_dados, prop_overlap, prop_sem_overlap) %>%
  gt() %>%
  tab_header(
    title = md("**RNA-Seq x STRs: Publication Summary**"),
    subtitle = "Per-gene counts and proportions of DBSCAN outliers and allele overlap"
  ) %>%
  cols_label(
    gene = "Gene",
    gse = "GSE Studies",
    n_strs_identified = "STRs in DEGs",
    n_strs_outliers = "With DBSCAN Outlier",
    prop_outliers = "% Outliers",
    n_overlap_sim = "Overlap (sim)",
    n_sem_dados = "No Data (sem dados)",
    prop_overlap = "% Overlap",
    prop_sem_overlap = "% No Overlap"
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_source_note(
    source_note = "Source: cross_DEGs_STRs.py output"
  ) %>%
  opt_stylize(style = 3)

out_gt_html <- file.path(out_dir, "rna_publication_table.html")
gtsave(pub_gt, out_gt_html)
cat(sprintf("Tabela gt salva em: %s\n", out_gt_html))

cat("\nConcluido.\n")
