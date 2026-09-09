#!/usr/bin/env Rscript
# plot_intervention_summary.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Gera visualizacoes para o cruzamento RNA-Seq x STRs POR INTERVENCAO:
#     1) Raincloud plot: distribuicao de allele2_est por intervencao,
#        colorido por group (case/control), com outliers DBSCAN.
#     2) Tabela de publicacao: por intervencao, contagens e proporcoes.
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog       STRs_analysis_dataset.tsv (per-sample x STR)
#   --intervention-strs intervention_strs.tsv (cruzação por intervenção)
#   --intervention-outliers intervention_outliers.tsv
#   --intervention-summary intervention_summary.tsv (resumo intervenção x gene)
#   --out-dir           Diretório de saída
#
# SAIDAS
#   intervention_raincloud.png
#   intervention_publication_table.tsv
#   intervention_publication_table.html
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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

path_str_catalog    <- parse_arg("--str-catalog")
path_intv_strs      <- parse_arg("--intervention-strs")
path_intv_outliers  <- parse_arg("--intervention-outliers")
path_intv_summary   <- parse_arg("--intervention-summary")
out_dir             <- parse_arg("--out-dir", ".")

missing_args <- c(
  "--str-catalog" = is.null(path_str_catalog),
  "--intervention-strs" = is.null(path_intv_strs),
  "--intervention-outliers" = is.null(path_intv_outliers),
  "--intervention-summary" = is.null(path_intv_summary)
)
if (any(missing_args)) {
  stop("Argumentos ausentes: ", paste(names(missing_args[missing_args]), collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("Carregando dados...\n")

str_cat <- fread(path_str_catalog, header = TRUE, sep = "\t")
cat(sprintf("  STR catalog: %d linhas\n", nrow(str_cat)))

intv_strs <- fread(path_intv_strs, header = TRUE, sep = "\t")
cat(sprintf("  intervention_strs.tsv: %d linhas\n", nrow(intv_strs)))

intv_out <- fread(path_intv_outliers, header = TRUE, sep = "\t")
cat(sprintf("  intervention_outliers.tsv: %d linhas\n", nrow(intv_out)))

intv_sum <- fread(path_intv_summary, header = TRUE, sep = "\t")
cat(sprintf("  intervention_summary.tsv: %d linhas\n", nrow(intv_sum)))

# ==========================================
# 2. Prepare data for raincloud
# ==========================================
cat("\nPreparando dados para raincloud...\n")

# Mapear gene_name -> intervention a partir de intv_strs
# Cada linha de intv_strs tem intervention, gene_name, STRs_ID, allele2_est, group
plot_data <- intv_strs[!is.na(group) & group != ""]
plot_data[, group := factor(group, levels = c("case", "control"))]
plot_data[, intervention := factor(intervention)]

cat(sprintf("  Intervenções: %s\n", paste(levels(plot_data$intervention), collapse = ", ")))
cat(sprintf("  Observações: %d (genes: %d)\n", nrow(plot_data), uniqueN(plot_data$gene_name)))

# ==========================================
# 3. Raincloud plot (scientific theme)
# ==========================================
cat("\nGerando raincloud plot por intervenção...\n")

group_colors <- c("case" = "#E41A1C", "control" = "#377EB8")

n_interv <- length(levels(plot_data$intervention))

p_rain <- ggplot(plot_data,
                 aes(x = allele2_est, y = intervention, fill = group)) +
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
    title = "Allele 2 length distributions by intervention",
    subtitle = "DBSCAN global outliers, colored by case/control status",
    x = "Allele 2 length (repeat units, log1p)",
    y = "Intervention",
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
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold",
                              margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, color = "grey40",
                                 margin = margin(b = 10)),
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

out_png <- file.path(out_dir, "intervention_raincloud.png")
raincloud_h <- max(5, n_interv * 0.8 + 3)
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
# 4. Publication table (per intervention)
# ==========================================
cat("\nGerando tabela de publicacao por intervenção...\n")

# --- 4.1 Outlier counts per intervention x group (sample-level) ---
outlier_by_intv <- intv_out[!is.na(group), .(
  n_outliers = .N,
  n_out_case = sum(group == "case"),
  n_out_control = sum(group == "control")
), by = intervention]

# --- 4.2 Outliers with no overlap (sample-level) ---
# Obter overlap_maior_alealo_grupos de intv_sum por gene
# Mapear gene -> overlap por intervention
overlap_map <- intv_sum[, .(intervention, gene, overlap_maior_alealo_grupos)]
intv_out_merged <- merge(intv_out, overlap_map, by = c("intervention", "gene"),
                         all.x = TRUE, allow.cartesian = TRUE)

out_no_overlap <- intv_out_merged[
  overlap_maior_alealo_grupos == "nao" & !is.na(group), .(
  n_out_no_overlap = .N,
  n_out_no_overlap_case = sum(group == "case"),
  n_out_no_overlap_control = sum(group == "control")
), by = intervention]

# --- 4.3 STR density per gene per intervention ---
strs_per_gene_intv <- intv_strs[, .(n_strs = uniqueN(STRs_ID)), by = .(intervention, gene_name)]
str_density_intv <- strs_per_gene_intv[, .(
  str_density_per_gene = round(sum(n_strs) / uniqueN(gene_name), 2)
), by = intervention]

# --- 4.4 Summary counts from intv_sum ---
summary_intv <- intv_sum[, .(
  n_genes = uniqueN(gene),
  n_strs_total = sum(n_strs_identified),
  n_overlap_sim = sum(overlap_maior_alealo_grupos == "sim", na.rm = TRUE),
  n_overlap_nao = sum(overlap_maior_alealo_grupos == "nao", na.rm = TRUE),
  n_sem_dados = sum(overlap_maior_alealo_grupos == "sem_dados", na.rm = TRUE)
), by = intervention]

# --- 4.5 Case/control counts from intv_strs ---
group_counts_intv <- intv_strs[!is.na(group), .(
  n_case = sum(group == "case"),
  n_control = sum(group == "control")
), by = intervention]

# --- 4.6 Median allele2_est per group per intervention ---
allele_medians_intv <- intv_strs[!is.na(group), .(
  median_allele_case = round(median(allele2_est[group == "case"], na.rm = TRUE), 2),
  median_allele_control = round(median(allele2_est[group == "control"], na.rm = TRUE), 2)
), by = intervention]

# --- 4.7 GSE info per intervention ---
gse_info <- unique(intv_strs[, .(intervention, gse)])

# --- 4.8 Merge all ---
pub_table <- merge(summary_intv, str_density_intv, by = "intervention", all.x = TRUE)
pub_table <- merge(pub_table, group_counts_intv, by = "intervention", all.x = TRUE)
pub_table <- merge(pub_table, allele_medians_intv, by = "intervention", all.x = TRUE)
pub_table <- merge(pub_table, outlier_by_intv, by = "intervention", all.x = TRUE)
pub_table <- merge(pub_table, out_no_overlap, by = "intervention", all.x = TRUE)
pub_table <- merge(pub_table, gse_info, by = "intervention", all.x = TRUE)

# Fill NA with 0
for (col in c("n_outliers", "n_out_case", "n_out_control",
              "n_out_no_overlap", "n_out_no_overlap_case", "n_out_no_overlap_control")) {
  pub_table[is.na(get(col)), (col) := 0]
}

# --- 4.9 Format columns: abs (pct%) ---
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

cat(sprintf("  Tabela final: %d intervenções\n", nrow(pub_table)))

# Save TSV
out_tsv <- file.path(out_dir, "intervention_publication_table.tsv")
pub_tsv <- pub_table[, .(intervention, gse, n_genes, n_strs_total, str_density_per_gene,
                          n_case_str, n_control_str,
                          median_allele_case, median_allele_control,
                          n_outliers_str, n_out_case_str, n_out_control_str,
                          n_overlap_str, n_out_no_overlap_str)]
cat("  Primeiras 3 linhas:\n")
print(head(as.data.frame(pub_tsv), 3))
fwrite(pub_tsv, out_tsv, sep = "\t")
cat(sprintf("Tabela salva em: %s (%d intervenções)\n", out_tsv, nrow(pub_tsv)))

# Save formatted gt table as HTML (publication ready)
pub_gt <- pub_tsv %>%
  gt() %>%
  tab_header(
    title = md("**Table X.** RNA-Seq STR analysis per intervention"),
    subtitle = "Summary of DEG STRs, DBSCAN outliers, and allele overlap by group"
  ) %>%
  tab_stubhead(label = "Intervention") %>%
  tab_spanner(
    label = "Study Information",
    columns = c(intervention, gse, n_genes, n_strs_total, str_density_per_gene)
  ) %>%
  tab_spanner(
    label = "Sample Distribution",
    columns = c(n_case_str, n_control_str)
  ) %>%
  tab_spanner(
    label = "Allele Statistics",
    columns = c(median_allele_case, median_allele_control)
  ) %>%
  tab_spanner(
    label = "DBSCAN Outliers",
    columns = c(n_outliers_str, n_out_case_str, n_out_control_str,
                n_overlap_str, n_out_no_overlap_str)
  ) %>%
  cols_label(
    intervention = "Intervention",
    gse = "GSE",
    n_genes = "Genes",
    n_strs_total = "STRs",
    str_density_per_gene = "Dens.",
    n_case_str = "Case",
    n_control_str = "Control",
    median_allele_case = "Median (case)",
    median_allele_control = "Median (ctrl)",
    n_outliers_str = "Total",
    n_out_case_str = "Case",
    n_out_control_str = "Control",
    n_overlap_str = "Overlap",
    n_out_no_overlap_str = "Sem overlap"
  ) %>%
  sub_missing(columns = everything(), missing_text = "-") %>%
  tab_style(
    style = list(
      cell_text(weight = "bold", size = px(11)),
      cell_borders(sides = "bottom", weight = px(1.5), color = "grey60")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(10)),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(size = px(10)),
    locations = cells_body()
  ) %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(10),
    heading.align = "left",
    column_labels.border.top.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.width = px(1),
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.width = px(1.5),
    table_body.border.bottom.color = "black",
    table_body.hlines.color = "grey90",
    table_body.hlines.width = px(0.5),
    table.border.left.width = px(0),
    table.border.right.width = px(0),
    data_row.padding = px(4),
    row_group.padding = px(6)
  ) %>%
  tab_source_note(
    source_note = "Values: absolute count (percentage%). Dens. = STRs per DEG gene."
  ) %>%
  tab_source_note(
    source_note = "Overlap = same allele range between case/control. Sem overlap = outlier observations without overlap."
  ) %>%
  tab_source_note(
    source_note = "Source: cross_intervention_STRs.py output"
  )

out_gt_html <- file.path(out_dir, "intervention_publication_table.html")
gtsave(pub_gt, out_gt_html)
cat(sprintf("Tabela gt HTML salva em: %s\n", out_gt_html))

cat("\nConcluido.\n")
