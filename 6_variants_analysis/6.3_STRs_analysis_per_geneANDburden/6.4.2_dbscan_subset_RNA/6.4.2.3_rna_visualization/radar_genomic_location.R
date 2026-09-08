#!/usr/bin/env Rscript
# radar_genomic_location.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Radar plots para visualizar a localizacao genomica de:
#     1) Outliers DBSCAN global (por regiao genômica)
#     2) Variantes sem sobreposicao do alelo maior entre grupos
#
#   Gera radar individual por GSE study + radar combinado.
#
# ENTRADAS
#   --rna-outliers      rna_outlier_genes.tsv (gene x STR x GSE, com region)
#   --rna-summary       rna_summary_by_study.tsv (gene x GSE, com overlap)
#   --out-dir           Diretorio de saida
#
# DEPENDENCIAS
#   R: data.table, ggplot2, fmsb, patchwork, scales
#   Instalar fmsb: install.packages("fmsb")
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

has_fmsb <- requireNamespace("fmsb", quietly = TRUE)
if (has_fmsb) {
  suppressPackageStartupMessages(library(fmsb))
  cat("Usando fmsb::radarchart()\n")
} else {
  cat("fmsb nao encontrado, usando ggplot2::coord_polar()\n")
}

# ==========================================
# Parse arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  return(default)
}

path_rna_outliers <- parse_arg("--rna-outliers")
path_rna_summary  <- parse_arg("--summary")
out_dir           <- parse_arg("--out-dir", ".")

missing <- c(
  "--rna-outliers" = is.null(path_rna_outliers),
  "--summary" = is.null(path_rna_summary)
)
if (any(missing)) {
  stop("Argumentos ausentes: ", paste(names(missing[missing]), collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("Carregando rna_outlier_genes.tsv...\n")
rna_outliers <- fread(path_rna_outliers, header = TRUE, sep = "\t")
cat(sprintf("  %d linhas (gene x STR x GSE)\n", nrow(rna_outliers)))

cat("Carregando rna_summary_by_study.tsv...\n")
rna_summary <- fread(path_rna_summary, header = TRUE, sep = "\t")
cat(sprintf("  %d linhas (gene x GSE)\n", nrow(rna_summary)))

# ==========================================
# 2. Prepare data: add overlap status at STR level
# ==========================================
cat("\nPreparando dados...\n")

# Merge overlap status from summary to outlier genes (by gse + gene)
str_data <- merge(
  rna_outliers,
  rna_summary[, .(gse, gene, overlap_maior_alealo_grupos)],
  by = c("gse", "gene"),
  all.x = TRUE
)

# Define region categories (standardized)
region_order <- c("promoter", "five_prime_utr", "three_prime_utr",
                  "exon", "intron", "non_coding_exons",
                  "intergenic", "others")
str_data[, region := ifelse(region %in% region_order, region, "others")]
str_data[, region := factor(region, levels = region_order)]

# Get all GSEs
all_gses <- sort(unique(str_data$gse))
cat(sprintf("  GSEs: %s\n", paste(all_gses, collapse = ", ")))

# ==========================================
# 3. Count STRs per region for each category
# ==========================================

# --- Outliers: STRs with n_outliers_dbscan_global >= 1
count_outliers_by_region <- function(dat) {
  outliers <- dat[!is.na(n_outliers_dbscan_global) & n_outliers_dbscan_global >= 1]
  counts <- outliers[, .N, by = region]
  # Ensure all regions present
  full <- data.table(region = region_order)
  counts <- merge(full, counts, by = "region", all.x = TRUE)
  counts[is.na(N), N := 0]
  return(counts)
}

# --- No overlap: STRs with overlap_maior_alealo_grupos == "nao"
count_nooverlap_by_region <- function(dat) {
  no_overlap <- dat[overlap_maior_alealo_grupos == "nao"]
  counts <- no_overlap[, .N, by = region]
  full <- data.table(region = region_order)
  counts <- merge(full, counts, by = "region", all.x = TRUE)
  counts[is.na(N), N := 0]
  return(counts)
}

# ==========================================
# 4. Radar plot function (fmsb or ggplot2 fallback)
# ==========================================
make_radar_fmsb <- function(counts, title, color, fill_alpha = 0.3) {
  # fmsb needs: max row, min row, then data row
  vals <- counts$N
  max_val <- max(vals, 1)
  radar_data <- data.frame(
    rbind(
      rep(max_val, length(vals)),   # max
      rep(0, length(vals)),          # min
      vals                           # values
    )
  )
  colnames(radar_data) <- counts$region
  rownames(radar_data) <- c("max", "min", "value")

  fmsb::radarchart(
    radar_data,
    axistype = 1,
    pcol = color, pfcol = scales::alpha(color, fill_alpha), plwd = 2, plty = 1,
    cglcol = "grey70", cglty = 1, cglwd = 0.8,
    axislabcol = "grey40",
    vlcex = 0.8,
    title = title
  )
}

make_radar_ggplot <- function(counts, title, color, fill_alpha = 0.3) {
  # ggplot2 coord_polar approach
  n <- nrow(counts)
  counts[, angle := 90 - 360 * (seq_len(n) - 0.5) / n]
  counts[, label := as.character(region)]

  ggplot(counts, aes(x = region, y = N, group = 1)) +
    geom_polygon(fill = scales::alpha(color, fill_alpha), color = color, linewidth = 0.8) +
    geom_point(color = color, size = 2) +
    coord_polar() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
}

# ==========================================
# 5. Generate radar plots per GSE + combined
# ==========================================
cat("\nGerando radar plots...\n")

# --- OUTLIERS ---
cat("  Outliers DBSCAN global...\n")

# Per GSE
radar_outlier_list <- list()
for (gse in all_gses) {
  dat_gse <- str_data[gse == gse]
  counts <- count_outliers_by_region(dat_gse)
  radar_outlier_list[[gse]] <- counts
  cat(sprintf("    %s: %d outliers em %d regioes\n",
              gse, sum(counts$N), sum(counts$N > 0)))
}

# Combined
counts_combined_outlier <- count_outliers_by_region(str_data)
cat(sprintf("    COMBINED: %d outliers em %d regioes\n",
            sum(counts_combined_outlier$N), sum(counts_combined_outlier$N > 0)))

# Plot per GSE
plot_list_outlier <- list()
for (gse in all_gses) {
  counts <- radar_outlier_list[[gse]]
  if (has_fmsb) {
    p <- ggplot() +
      annotation_custom(
        grob = grid::recordGrob(make_radar_fmsb(counts, gse, "#E41A1C")),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      theme_void()
  } else {
    p <- make_radar_ggplot(counts, gse, "#E41A1C")
  }
  plot_list_outlier[[gse]] <- p
}

# Combined plot
if (has_fmsb) {
  p_comb_outlier <- ggplot() +
    annotation_custom(
      grob = grid::recordGrob(make_radar_fmsb(counts_combined_outlier, "COMBINED", "#E41A1C")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    theme_void()
} else {
  p_comb_outlier <- make_radar_ggplot(counts_combined_outlier, "COMBINED", "#E41A1C")
}
plot_list_outlier[["COMBINED"]] <- p_comb_outlier

# Arrange and save
combined_outlier <- wrap_plots(plot_list_outlier, ncol = 2)
outlier_png <- file.path(out_dir, "radar_outliers_by_study.png")
ggsave(outlier_png, combined_outlier, width = 14, height = 10, dpi = 300, bg = "white")
cat(sprintf("  Salvo: %s\n", outlier_png))

# --- NO OVERLAP ---
cat("  Variantes sem sobreposicao...\n")

radar_nooverlap_list <- list()
for (gse in all_gses) {
  dat_gse <- str_data[gse == gse]
  counts <- count_nooverlap_by_region(dat_gse)
  radar_nooverlap_list[[gse]] <- counts
  cat(sprintf("    %s: %d sem sobreposicao em %d regioes\n",
              gse, sum(counts$N), sum(counts$N > 0)))
}

counts_combined_nooverlap <- count_nooverlap_by_region(str_data)
cat(sprintf("    COMBINED: %d sem sobreposicao em %d regioes\n",
            sum(counts_combined_nooverlap$N), sum(counts_combined_nooverlap$N > 0)))

# Plot per GSE
plot_list_nooverlap <- list()
for (gse in all_gses) {
  counts <- radar_nooverlap_list[[gse]]
  if (has_fmsb) {
    p <- ggplot() +
      annotation_custom(
        grob = grid::recordGrob(make_radar_fmsb(counts, gse, "#FFBF00")),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      theme_void()
  } else {
    p <- make_radar_ggplot(counts, gse, "#FFBF00")
  }
  plot_list_nooverlap[[gse]] <- p
}

# Combined
if (has_fmsb) {
  p_comb_nooverlap <- ggplot() +
    annotation_custom(
      grob = grid::recordGrob(make_radar_fmsb(counts_combined_nooverlap, "COMBINED", "#FFBF00")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    theme_void()
} else {
  p_comb_nooverlap <- make_radar_ggplot(counts_combined_nooverlap, "COMBINED", "#FFBF00")
}
plot_list_nooverlap[["COMBINED"]] <- p_comb_nooverlap

combined_nooverlap <- wrap_plots(plot_list_nooverlap, ncol = 2)
nooverlap_png <- file.path(out_dir, "radar_no_overlap_by_study.png")
ggsave(nooverlap_png, combined_nooverlap, width = 14, height = 10, dpi = 300, bg = "white")
cat(sprintf("  Salvo: %s\n", nooverlap_png))

# ==========================================
# 6. Summary table
# ==========================================
cat("\nGerando tabela resumo...\n")

summary_rows <- list()

# Outliers per region per GSE
for (gse in all_gses) {
  counts <- radar_outlier_list[[gse]]
  for (i in seq_len(nrow(counts))) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      gse = gse, region = as.character(counts$region[i]),
      category = "outliers", count = counts$N[i]
    )
  }
}
# Combined outliers
for (i in seq_len(nrow(counts_combined_outlier))) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    gse = "COMBINED", region = as.character(counts_combined_outlier$region[i]),
    category = "outliers", count = counts_combined_outlier$N[i]
  )
}
# No overlap per GSE
for (gse in all_gses) {
  counts <- radar_nooverlap_list[[gse]]
  for (i in seq_len(nrow(counts))) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      gse = gse, region = as.character(counts$region[i]),
      category = "no_overlap", count = counts$N[i]
    )
  }
}
# Combined no overlap
for (i in seq_len(nrow(counts_combined_nooverlap))) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    gse = "COMBINED", region = as.character(counts_combined_nooverlap$region[i]),
    category = "no_overlap", count = counts_combined_nooverlap$N[i]
  )
}

summary_table <- rbindlist(summary_rows)
out_tsv <- file.path(out_dir, "radar_genomic_summary.tsv")
fwrite(summary_table, out_tsv, sep = "\t")
cat(sprintf("  Tabela: %s (%d linhas)\n", out_tsv, nrow(summary_table)))

# ==========================================
# 7. Final summary
# ==========================================
cat("\n", strrep("=", 55), "\n")
cat(" RESUMO: Radar plots — localizacao genômica\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  GSEs analisados: %d (%s)\n", length(all_gses),
            paste(all_gses, collapse = ", ")))
cat(sprintf("  Regioes genômicas: %s\n", paste(region_order, collapse = ", ")))
cat(sprintf("  Total STRs analisados: %d\n", nrow(str_data)))
cat(sprintf("  Outliers DBSCAN global: %d\n",
            sum(!is.na(str_data$n_outliers_dbscan_global) &
                    str_data$n_outliers_dbscan_global >= 1)))
cat(sprintf("  Sem sobreposicao: %d\n",
            sum(str_data$overlap_maior_alealo_grupos == "nao", na.rm = TRUE)))
cat(sprintf("  Library usada: %s\n", ifelse(has_fmsb, "fmsb", "ggplot2::coord_polar")))
cat(strrep("=", 55), "\n")
cat("Concluido.\n")
