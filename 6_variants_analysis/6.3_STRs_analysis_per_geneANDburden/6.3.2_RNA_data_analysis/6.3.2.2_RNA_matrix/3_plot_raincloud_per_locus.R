#!/usr/bin/env Rscript
# plot_raincloud_per_locus.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Raincloud POR LOCUS (STR) com apenas outliers DBSCAN.
#   Para cada intervenção (ou uma única, se --intervention) gera, para
#   cada GSE (acesso do estudo) e para a visão agregada (ALL):
#     1) Raincloud por locus, salvo em <out-dir>/<GSE>/<intervention>_*.png
#     2) CSV por paciente, salvo em <out-dir>/<GSE>/<intervention>_patients.csv
#
# ENTRADAS
#   --intervention-outliers intervention_outliers.tsv (outliers DBSCAN)
#   --intervention          <nome> | "ALL" (default: ALL)
#   --out-dir               Diretorio de saida
#
# SAIDAS
#   <out-dir>/ALL/<intervention>_patients.csv
#   <out-dir>/ALL/<intervention>_raincloud_per_locus_pXX.png   (facet todos GSEs)
#   <out-dir>/<GSE>/<intervention>_patients.csv
#   <out-dir>/<GSE>/<intervention>_raincloud_per_locus_pXX.png
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrain)
  library(scales)
})

options(ragg.max_dim = 200000)

# ==========================================
# Parse arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  return(default)
}

path_intv_outliers <- parse_arg("--intervention-outliers")
intervention_sel   <- parse_arg("--intervention", "ALL")
out_dir            <- parse_arg("--out-dir", ".")

if (is.null(path_intv_outliers)) {
  stop("Argumento ausente: --intervention-outliers")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_colors <- c("case" = "#E41A1C", "control" = "#377EB8")
loci_per_page <- 20

# ==========================================
# 1. Load data
# ==========================================
cat("Carregando intervention_outliers.tsv...\n")
intv_out <- fread(path_intv_outliers, header = TRUE, sep = "\t")
cat(sprintf("  Linhas: %d\n", nrow(intv_out)))

for (req in c("STRs_ID", "intervention", "group", "gene_name", "gse"))
  if (!req %in% names(intv_out)) stop("Coluna ausente no input: ", req)

intv_out[, group := tolower(trimws(as.character(group)))]
intv_out <- intv_out[!is.na(group) & group != ""]
intv_out[, group := factor(group, levels = c("case", "control"))]

intv_out[, locus_label := ifelse(is.na(gene_name) | gene_name == "",
                                 STRs_ID,
                                 sprintf("%s | %s", gene_name, STRs_ID))]

available <- sort(unique(intv_out$intervention))
cat(sprintf("  Intervenções disponíveis: %s\n",
            paste(available, collapse = ", ")))

if (toupper(intervention_sel) == "ALL") {
  interventions <- available
} else {
  if (!intervention_sel %in% available) {
    stop("Intervenção '", intervention_sel, "' não encontrada. Disponíveis: ",
         paste(available, collapse = ", "))
  }
  interventions <- intervention_sel
}

# ==========================================
# Helper: gera rainclouds de um dataset (paginação por locus)
# ==========================================
# gse_tag: "ALL" (facet por gse) ou nome do GSE. Define a pasta de destino.
make_raincloud <- function(dat, intv, gse_tag, title, subtitle) {
  dest_dir <- file.path(out_dir, gsub("[^A-Za-z0-9._-]", "_", gse_tag))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  loci <- unique(dat$locus_label)
  n_loci <- length(loci)
  cat(sprintf("    [%s] %d loci -> %d páginas\n",
              gse_tag, n_loci, ceiling(n_loci / loci_per_page)))

  pages <- split(loci, ceiling(seq_along(loci) / loci_per_page))

  for (pg_i in seq_along(pages)) {
    page_loci <- pages[[pg_i]]
    pdat <- dat[locus_label %in% page_loci]
    pdat[, locus_label := factor(locus_label,
                                 levels = rev(sort(unique(locus_label))))]

    p <- ggplot(pdat,
                aes(x = allele2_est, y = locus_label, fill = group)) +
      geom_rain(
        alpha = 0.5,
        point.args = list(size = 1.2, alpha = 0.6),
        boxplot.args = list(outlier.shape = NA, width = 0.2)
      ) +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_x_continuous(
        trans = "log1p",
        breaks = c(0, 1, 5, 10, 50, 100, 500, 1000),
        labels = c("0", "1", "5", "10", "50", "100", "500", "1k"),
        expand = c(0.01, 0)
      ) +
      labs(
        title = title,
        subtitle = subtitle,
        x = "Allele length (repeat units, log1p)",
        y = "Locus (gene | STRs_ID)",
        caption = "Only DBSCAN global outliers shown"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 7, color = "grey20"),
        axis.text.x = element_text(size = 9, color = "grey20"),
        plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5),
        legend.position = "bottom",
        panel.spacing = unit(1.2, "lines"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    # ATENCAO: facet por gse apenas na visao agregada (ALL)
    if (toupper(gse_tag) == "ALL") {
      p <- p + facet_wrap(~ gse, scales = "free_y", ncol = 1) +
        theme(strip.text = element_text(size = 11, face = "bold"))
    }

    out_file <- file.path(dest_dir,
      sprintf("%s_raincloud_per_locus_p%02d.png",
              gsub("[^A-Za-z0-9._-]", "_", intv), pg_i))
    ggsave(
      filename = out_file,
      plot = p,
      width = 12,
      height = max(6, length(page_loci) * 0.35 + 3),
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    cat(sprintf("    Salvo: %s\n", out_file))
  }
}

# ==========================================
# Helper: CSV por paciente numa pasta de GSE
# ==========================================
write_patients_csv <- function(dat, intv, gse_tag) {
  dest_dir <- file.path(out_dir, gsub("[^A-Za-z0-9._-]", "_", gse_tag))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  csv_out <- file.path(dest_dir,
                       sprintf("%s_patients.csv",
                               gsub("[^A-Za-z0-9._-]", "_", intv)))
  csv_cols <- c("STRs_ID", "locus_label", "gene_name", "chrom", "start",
                "end", "repeat_unit", "region", "allele1_est", "allele2_est",
                "depth", "group", "gse", "dataset",
                "logFC", "FDR", "P.Value", "Direction",
                "n_outliers_dbscan_global", "outlier_samples_dbscan_global",
                "n_clusters_dbscan_global", "noise_ratio_dbscan_global")
  present_cols <- intersect(csv_cols, names(dat))
  fwrite(dat[, ..present_cols], csv_out)
  cat(sprintf("    CSV por paciente: %s\n", csv_out))
}

# ==========================================
# 2. Loop por intervenção
# ==========================================
for (intv in interventions) {
  cat(sprintf("\n=== Intervenção: %s ===\n", intv))

  dat <- intv_out[intervention == intv]
  if (nrow(dat) == 0) {
    cat("  Sem observações, pulando.\n")
    next
  }
  cat(sprintf("  Observações: %d | Loci: %d | GSEs: %s\n",
              nrow(dat), uniqueN(dat$locus_label),
              paste(sort(unique(dat$gse)), collapse = ", ")))

  # --- 2.1 Visão agregada (ALL): CSV + raincloud facet por gse ---
  cat("  Gerando visão agregada (todos GSEs)...\n")
  write_patients_csv(dat, intv, "ALL")
  make_raincloud(
    dat = dat, intv = intv, gse_tag = "ALL",
    title = sprintf("DBSCAN outliers per locus - %s", gsub("_", " ", intv)),
    subtitle = "All GSEs | alleles by case/control"
  )

  # --- 2.2 Por GSE (acesso do estudo): CSV + raincloud ---
  for (gse_i in sort(unique(dat$gse))) {
    dat_gse <- dat[gse == gse_i]
    cat(sprintf("  Gerando por GSE: %s\n", gse_i))
    write_patients_csv(dat_gse, intv, gse_i)
    make_raincloud(
      dat = dat_gse, intv = intv, gse_tag = gse_i,
      title = sprintf("DBSCAN outliers per locus - %s (%s)",
                      gsub("_", " ", intv), gse_i),
      subtitle = sprintf("%s | alleles by case/control", gse_i)
    )
  }
}

cat("\nConcluido.\n")