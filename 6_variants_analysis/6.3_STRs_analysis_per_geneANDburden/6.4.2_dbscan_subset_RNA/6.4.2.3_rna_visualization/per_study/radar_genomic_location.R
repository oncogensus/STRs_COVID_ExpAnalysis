#!/usr/bin/env Rscript
# radar_genomic_location.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Radar plots (coord_polar) para visualizar a localizacao genomica de:
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
#   R: data.table, ggplot2, patchwork, scales
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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

str_data <- merge(
  rna_outliers,
  rna_summary[, .(gse, gene, overlap_maior_alealo_grupos)],
  by = c("gse", "gene"),
  all.x = TRUE
)

region_order <- c("promoter", "five_prime_utr", "three_prime_utr",
                  "exon", "intron", "non_coding_exons",
                  "intergenic", "others")
str_data[, region := ifelse(region %in% region_order, region, "others")]
str_data[, region := factor(region, levels = region_order)]

all_gses <- sort(unique(str_data$gse))
cat(sprintf("  GSEs: %s\n", paste(all_gses, collapse = ", ")))

# ==========================================
# 3. Count STRs per region
# ==========================================
count_by_region <- function(dat, filter_fn) {
  filtered <- dat[filter_fn(dat)]
  counts <- filtered[, .N, by = region]
  full <- data.table(region = region_order)
  counts <- merge(full, counts, by = "region", all.x = TRUE)
  counts[is.na(N), N := 0]
  return(counts)
}

# ==========================================
# 4. Radar plot function (ggplot2 coord_polar)
# ==========================================
make_radar <- function(counts, title, color, fill_alpha = 0.35) {
  n <- nrow(counts)
  counts[, label := as.character(region)]
  counts[, angle := 90 - 360 * (seq_len(n) - 0.5) / n]

  max_val <- max(counts$N, 1)

  ggplot(counts, aes(x = region, y = N, group = 1)) +
    geom_polygon(
      fill = scales::alpha(color, fill_alpha),
      color = color, linewidth = 0.9, lineend = "round"
    ) +
    geom_point(color = color, size = 2.5) +
    geom_text(
      aes(label = ifelse(N > 0, as.character(N), "")),
      vjust = -1.2, size = 3, color = "grey30"
    ) +
    scale_y_continuous(limits = c(0, max_val * 1.2), expand = c(0, 0)) +
    coord_polar() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7, color = "grey30"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5,
                                margin = margin(b = 5))
    )
}

# ==========================================
# 5. Generate radar plots
# ==========================================
cat("\nGerando radar plots...\n")

# --- OUTLIERS ---
cat("  Outliers DBSCAN global...\n")

radar_outlier_list <- list()
for (gse in all_gses) {
  counts <- count_by_region(str_data[gse == gse],
                            function(d) !is.na(d$n_outliers_dbscan_global) & d$n_outliers_dbscan_global >= 1)
  radar_outlier_list[[gse]] <- counts
  cat(sprintf("    %s: %d outliers em %d regioes\n",
              gse, sum(counts$N), sum(counts$N > 0)))
}

counts_comb_out <- count_by_region(str_data,
                                   function(d) !is.na(d$n_outliers_dbscan_global) & d$n_outliers_dbscan_global >= 1)
cat(sprintf("    COMBINED: %d outliers em %d regioes\n",
            sum(counts_comb_out$N), sum(counts_comb_out$N > 0)))

# Build plots
plot_list_out <- list()
for (gse in all_gses) {
  plot_list_out[[gse]] <- make_radar(radar_outlier_list[[gse]], gse, "#E41A1C")
}
plot_list_out[["COMBINED"]] <- make_radar(counts_comb_out, "COMBINED", "#C62828")

combined_outlier <- wrap_plots(plot_list_out, ncol = 2)
outlier_png <- file.path(out_dir, "radar_outliers_by_study.png")
ggsave(outlier_png, combined_outlier, width = 14, height = 10, dpi = 300, bg = "white")
cat(sprintf("  Salvo: %s\n", outlier_png))

# --- NO OVERLAP ---
cat("  Variantes sem sobreposicao...\n")

radar_nooverlap_list <- list()
for (gse in all_gses) {
  counts <- count_by_region(str_data[gse == gse],
                            function(d) d$overlap_maior_alealo_grupos == "nao")
  radar_nooverlap_list[[gse]] <- counts
  cat(sprintf("    %s: %d sem sobreposicao em %d regioes\n",
              gse, sum(counts$N), sum(counts$N > 0)))
}

counts_comb_no <- count_by_region(str_data,
                                  function(d) d$overlap_maior_alealo_grupos == "nao")
cat(sprintf("    COMBINED: %d sem sobreposicao em %d regioes\n",
            sum(counts_comb_no$N), sum(counts_comb_no$N > 0)))

plot_list_no <- list()
for (gse in all_gses) {
  plot_list_no[[gse]] <- make_radar(radar_nooverlap_list[[gse]], gse, "#FFBF00")
}
plot_list_no[["COMBINED"]] <- make_radar(counts_comb_no, "COMBINED", "#F57F17")

combined_nooverlap <- wrap_plots(plot_list_no, ncol = 2)
nooverlap_png <- file.path(out_dir, "radar_no_overlap_by_study.png")
ggsave(nooverlap_png, combined_nooverlap, width = 14, height = 10, dpi = 300, bg = "white")
cat(sprintf("  Salvo: %s\n", nooverlap_png))

# ==========================================
# 6. Summary table
# ==========================================
cat("\nGerando tabela resumo...\n")

summary_rows <- list()
for (gse in c(all_gses, "COMBINED")) {
  dat_src <- if (gse == "COMBINED") str_data else str_data[gse == gse]
  # Outliers
  co <- if (gse == "COMBINED") counts_comb_out else radar_outlier_list[[gse]]
  for (i in seq_len(nrow(co))) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      gse = gse, region = as.character(co$region[i]),
      category = "outliers", count = co$N[i]
    )
  }
  # No overlap
  cn <- if (gse == "COMBINED") counts_comb_no else radar_nooverlap_list[[gse]]
  for (i in seq_len(nrow(cn))) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      gse = gse, region = as.character(cn$region[i]),
      category = "no_overlap", count = cn$N[i]
    )
  }
}

summary_table <- rbindlist(summary_rows)
out_tsv <- file.path(out_dir, "radar_genomic_summary.tsv")
fwrite(summary_table, out_tsv, sep = "\t")
cat(sprintf("  Tabela: %s (%d linhas)\n", out_tsv, nrow(summary_table)))

# ==========================================
# 7. Summary
# ==========================================
cat("\n", strrep("=", 55), "\n")
cat(" RESUMO: Radar plots — localizacao genômica\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  GSEs analisados: %d (%s)\n", length(all_gses),
            paste(all_gses, collapse = ", ")))
cat(sprintf("  Regioes: %s\n", paste(region_order, collapse = ", ")))
cat(sprintf("  Total STRs: %d\n", nrow(str_data)))
cat(sprintf("  Outliers DBSCAN global: %d\n",
            sum(!is.na(str_data$n_outliers_dbscan_global) &
                    str_data$n_outliers_dbscan_global >= 1)))
cat(sprintf("  Sem sobreposicao: %d\n",
            sum(str_data$overlap_maior_alealo_grupos == "nao", na.rm = TRUE)))
cat(strrep("=", 55), "\n")
cat("Concluido.\n")
