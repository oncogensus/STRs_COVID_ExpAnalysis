#!/usr/bin/env Rscript
# plot_intervention_summary.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Gera visualizacoes para o cruzamento RNA-Seq x STRs POR INTERVENCAO:
#     1) 2D hexbin density: allele2_est vs depth por GSE e group.
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
#   intervention_density.png
#   intervention_publication_table.tsv
#   intervention_publication_table.html
#
# ESTILO
#   Density plot em hexbin com paleta Spectral (RColorBrewer), inspirado
#   na estetica classica do pacote {hexbin}:
#     bin <- hexbin(x, y, xbins = 40)
#     my_colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
#     plot(bin, colramp = my_colors, legend = FALSE)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gt)
  library(scales)
  library(hexbin)
  library(RColorBrewer)
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
# 2. Prepare data for density plot
# ==========================================
cat("\nPreparando dados para density plot...\n")

plot_data <- intv_out[!is.na(group) & group != ""]
plot_data[, group := factor(group, levels = c("case", "control"))]

cat(sprintf("  Intervenções: %s\n", paste(unique(plot_data$intervention), collapse = ", ")))
cat(sprintf("  Observações: %d (GSEs: %d)\n", nrow(plot_data), uniqueN(plot_data$gse)))

# ==========================================
# 3. 2D HEXBIN density: allele2_est vs depth
# ==========================================
cat("\nGerando density plot (hexbin)...\n")

n_gse <- uniqueN(plot_data$gse)

# ---- Paleta Spectral revertida, como no exemplo do {hexbin} -----------
#      rev(brewer.pal(11, "Spectral")) => azul (baixo) -> vermelho (alto)
hex_colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)

# ---- Numero de bins adaptativo: 40 para datasets grandes ---------------
#      (o exemplo original usa xbins = 40; em facetas com muitos pontos
#      um pouco menos de bins ajuda a manter o grafico legivel)
hex_bins <- 40

p_density <- ggplot(plot_data, aes(x = allele2_est, y = depth)) +
  # ---- Hexbin em vez de stat_density_2d ----
  geom_hex(bins = hex_bins, colour = NA) +
  # ---- Paleta Spectral revertida ----
  scale_fill_gradientn(
    colors = hex_colors,
    name   = "Count",
    guide  = guide_colorbar(
      barwidth  = unit(12, "mm"),
      barheight = unit(80, "mm"),
      title.position = "top"
    )
  ) +
  facet_wrap(~ gse + group, scales = "free", ncol = 2) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title    = "2D Density: Allele size vs Coverage",
    subtitle = "DBSCAN global outliers, faceted by GSE and group",
    x        = "Allele 2 length (repeat units)",
    y        = "Coverage (depth)",
    caption  = "DBSCAN global outliers only"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    # Grid removido dentro das facetas para nao competir com os hexbinos
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "grey20"),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold",
                              margin = margin(b = 15)),
    plot.subtitle = element_text(size = 11, color = "grey40",
                                 margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    panel.spacing = unit(1.2, "lines"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

out_png <- file.path(out_dir, "intervention_density.png")
density_h <- max(5, n_gse * 3 + 2)
ggsave(
  filename = out_png,
  plot = p_density,
  width = 10,
  height = density_h,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
cat(sprintf("Density plot salvo em: %s\n", out_png))

# ==========================================
# 4. Publication table (per intervention)
# ==========================================
cat("\nGerando tabela de publicacao por intervenção...\n")

# --- 4.1 Outlier counts per intervention x group (sample-level) ---
setnames(intv_out, "gene_name", "gene")

outlier_by_intv <- intv_out[!is.na(group), .(
  n_outliers = .N,
  n_out_case = sum(group == "case"),
  n_out_control = sum(group == "control")
), by = intervention]

# --- 4.2 Outliers with no overlap (sample-level) ---
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