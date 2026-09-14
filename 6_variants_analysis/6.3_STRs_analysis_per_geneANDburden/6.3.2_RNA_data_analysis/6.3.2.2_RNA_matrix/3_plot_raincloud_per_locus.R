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
#
# ESTILO
#   Raincloud em ggplot2 puro (sem {ggrain}), seguindo a abordagem de
#   Cedric Scherer:
#   https://www.cedricscherer.com/2021/06/06/
#   visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/
#
#   GeomFlatViolin adaptado para ggplot2 >= 3.4 (xmin/xmax manuais).
#   Boxplot, pontos e violino sao deslocados com position_nudge para
#   evitar sobreposicao dentro de cada locus.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

options(ragg.max_dim = 200000)

# ==========================================
# geom_flat_violin — meia-violino para raincloud
# Definida manualmente (nao existe nativamente no ggplot2).
# Adaptada do post de Cedric Scherer, com fix para ggplot2 >= 3.4.
# ==========================================
"%||%" <- function(a, b) if (!is.null(a)) a else b

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomFlatViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, ...)
  )
}

GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,

  setup_data = function(data, params) {
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)

    # --- FIX ggplot2 >= 3.4 ---------------------------------------------
    # xmin/xmax deixaram de ser adicionados automaticamente pela escala
    # discreta. Calculamos manualmente para o draw_group funcionar.
    # --------------------------------------------------------------------
    data$x    <- as.numeric(data$x)
    data$xmin <- data$x - data$width / 2
    data$xmax <- data$x + data$width / 2

    data
  },

  draw_group = function(data, panel_params, coord) {
    data <- transform(data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x))
    newdata <- rbind(
      transform(data, x = xminv)[order(data$y), ],
      transform(data, x = xmaxv)[order(data$y, decreasing = TRUE), ]
    )
    newdata <- rbind(newdata, newdata[1, ])
    newdata$group <- 1
    ggplot2:::ggname("geom_flat_violin",
                     GeomPolygon$draw_panel(newdata, panel_params, coord))
  },

  draw_key = draw_key_polygon,
  default_aes = aes(weight = 1, colour = "grey20", fill = "white",
                    linewidth = 0.5, linetype = "solid", alpha = NA),
  required_aes = c("x", "y")
)

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

group_colors  <- c("case" = "#E41A1C", "control" = "#377EB8")
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

intv_out[, str_variant := sub("^[^:]+:[^:]+:(.+)$", "\\1", STRs_ID)]
intv_out[, locus_label := ifelse(is.na(gene_name) | gene_name == "",
                                 str_variant,
                                 sprintf("%s (%s)", gene_name, str_variant))]

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
#
# LAYOUT (dentro de um locus, apos coord_flip):
#   x-0.36 .. x-0.24  -> boxplot  (nudge -0.30, width 0.12)
#   x-0.03 .. x+0.03  -> pontos   (jitter 0.03)
#   x+0.05 .. x+0.55  -> violino  (nudge +0.30, width 0.50)
# Assim os tres elementos nao se sobrepoem.
# ==========================================
make_raincloud <- function(dat, intv, gse_tag, title, subtitle) {
  dest_dir <- file.path(out_dir, gsub("[^A-Za-z0-9._-]", "_", gse_tag))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  loci   <- unique(dat$locus_label)
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
                aes(x = locus_label, y = allele2_est,
                    fill = group, colour = group)) +
      # ---- 1) Meia-violino (a "nuvem"), a direita do locus ----
      geom_flat_violin(
        position = position_nudge(x = 0.30),
        trim     = FALSE,
        alpha    = 0.5,
        colour   = NA,
        width    = 0.50
      ) +
      # ---- 2) Boxplot (a "chuva" inferior), a esquerda do locus ----
      geom_boxplot(
        width         = 0.12,
        outlier.shape = NA,
        alpha         = 0.6,
        colour        = "black",
        position      = position_nudge(x = -0.30),
        show.legend   = FALSE
      ) +
      # ---- 3) Pontos jittered (a "chuva" superior), no centro do locus ----
      geom_point(
        size        = 3,
        alpha       = 0.5,
        position    = position_jitter(width = 0.03, seed = 321),
        show.legend = FALSE
      ) +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_colour_manual(values = group_colors, guide = "none") +
      # coord_flip deixa o raincloud horizontal (locus no eixo Y)
      coord_flip() +
      scale_y_continuous(
        trans  = "log1p",
        breaks = c(0, 1, 5, 10, 50, 100, 500, 1000),
        labels = c("0", "1", "5", "10", "50", "100", "500", "1k"),
        expand = c(0.01, 0)
      ) +
      labs(
        title    = title,
        subtitle = subtitle,
        x        = NULL,
        y        = "Allele length (repeat units, log1p)",
      ) +
      theme_classic(base_size = 11) +
      theme(
        # Eixo categórico (locus) limpo
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.y      = element_text(size = 7, color = "grey20"),
        # Eixo contínuo (allele length)
        axis.text.x      = element_text(size = 9, color = "grey20"),
        axis.title.x     = element_text(size = 11, face = "bold"),
        # Títulos
        plot.title       = element_text(size = 13, hjust = 0.5, face = "bold"),
        plot.subtitle    = element_text(size = 10, color = "grey40", hjust = 0.5),
        # Strip (facet) sutil
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(size = 11, face = "bold"),
        # Legenda embaixo
        legend.position  = "bottom",
        legend.title     = element_text(size = 10),
        # Espaçamento entre facets
        panel.spacing    = unit(1.2, "lines"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)
      )

    # ATENCAO: facet por gse apenas na visao agregada (ALL)
    if (toupper(gse_tag) == "ALL") {
      p <- p + facet_wrap(~ gse, scales = "free_y", ncol = 1)
    }

    out_file <- file.path(dest_dir,
      sprintf("%s_raincloud_per_locus_p%02d.png",
              gsub("[^A-Za-z0-9._-]", "_", intv), pg_i))
    ggsave(
      filename  = out_file,
      plot      = p,
      width     = 12,
      height    = max(6, length(page_loci) * 0.35 + 3),
      dpi       = 600,
      bg        = "white",
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