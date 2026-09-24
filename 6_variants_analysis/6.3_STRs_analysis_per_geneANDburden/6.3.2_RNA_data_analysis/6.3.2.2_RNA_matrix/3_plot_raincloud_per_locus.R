#!/usr/bin/env Rscript
# plot_raincloud_per_locus.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Raincloud POR LOCUS (STR) com apenas outliers DBSCAN.
#   Para cada intervenção (ou uma única, se --intervention) gera, para
#   cada GSE (acesso do estudo):
#     1) Raincloud por locus, salvo em <out-dir>/<GSE>/<intervention>_*.png
#     2) CSV por paciente, salvo em <out-dir>/<GSE>/<intervention>_patients.csv
#
# ENTRADAS
#   --intervention-outliers intervention_outliers.tsv (outliers DBSCAN)
#   --intervention          <nome> | "ALL" (default: ALL)
#   --out-dir               Diretorio de saida
#
# SAIDAS
#   <out-dir>/<GSE>/<intervention>_patients.csv
#   <out-dir>/<GSE>/<intervention>_raincloud_per_locus_pXX.png
#
# ESTILO
#   Raincloud em ggplot2 puro. Repeat units em escala LINEAR (sem log1p).
#   Sem subtitle. Exceção: COVID_vs_CONTROL||GSE183533 (split + box/dots
#   separados + caption com o critério de estratificação).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

options(ragg.max_dim = 200000)

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
# Overrides de estetica por (intervencao, GSE)
# ==========================================
style_overrides <- list(
  "COVID_vs_CONTROL||GSE183533" = list(
    box_dots_separate = TRUE,
    split_by_allele   = TRUE,
    show_violin       = FALSE,
    box_nudge         = -0.22,
    box_width         = 0.20,
    dots_nudge        = 0.22,
    point_jitter      = 0.05,
    point_size        = 2.5,
    height_factor     = 0.35,
    height_offset     = 3
  )
)

default_style <- list(
  box_dots_separate = FALSE,
  split_by_allele   = FALSE,
  show_violin       = TRUE,
  violin_nudge      = 0.30,
  violin_width      = 0.50,
  box_nudge         = -0.30,
  box_width         = 0.12,
  dots_nudge        = 0.00,
  point_size        = 3,
  point_jitter      = 0.03,
  height_factor     = 0.35,
  height_offset     = 3
)

get_style <- function(intv, gse_tag) {
  key_exact <- paste0(intv, "||", gse_tag)
  key_any   <- paste0(intv, "||*")
  if (!is.null(style_overrides[[key_exact]])) return(style_overrides[[key_exact]])
  if (!is.null(style_overrides[[key_any]]))   return(style_overrides[[key_any]])
  default_style
}

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
make_raincloud <- function(dat, intv, gse_tag, title, subtitle) {
  dest_dir <- file.path(out_dir, gsub("[^A-Za-z0-9._-]", "_", gse_tag))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  loci   <- unique(dat$locus_label)
  n_loci <- length(loci)
  cat(sprintf("    [%s] %d loci -> %d páginas\n",
              gse_tag, n_loci, ceiling(n_loci / loci_per_page)))

  pages <- split(loci, ceiling(seq_along(loci) / loci_per_page))

  st <- get_style(intv, gse_tag)
  cat(sprintf("    [%s] estilo: sep=%s split=%s violin=%s\n",
              gse_tag, st$box_dots_separate, st$split_by_allele,
              st$show_violin))

  for (pg_i in seq_along(pages)) {
    page_loci <- pages[[pg_i]]
    pdat <- dat[locus_label %in% page_loci]
    pdat[, locus_label := factor(locus_label,
                                 levels = rev(sort(unique(locus_label))))]

    # ---- Split por tamanho do alelo (opcional) ----
    use_facet     <- FALSE
    split_caption <- NULL

    if (isTRUE(st$split_by_allele) && n_loci > 2) {
      loci_max <- pdat[, .(max_al = max(allele2_est, na.rm = TRUE)),
                       by = locus_label]
      thr <- median(loci_max$max_al, na.rm = TRUE)

      loci_max[, size_group := ifelse(max_al <= thr, "Short STRs", "Long STRs")]
      pdat <- merge(pdat, loci_max[, .(locus_label, size_group)],
                    by = "locus_label", all.x = TRUE)
      pdat[, size_group := factor(size_group,
                                   levels = c("Short STRs", "Long STRs"))]
      setorder(loci_max, max_al)
      pdat[, locus_label := factor(locus_label, levels = loci_max$locus_label)]
      use_facet <- TRUE

      split_caption <- sprintf(
        paste0("Loci were stratified into Short and Long STRs by the median ",
               "of the per-locus maximum allele length ",
               "(threshold = %d repeat units)."),
        round(thr)
      )

      cat(sprintf("    [%s] split: %d short + %d long (thr=%d)\n",
                  gse_tag,
                  sum(loci_max$max_al <= thr), sum(loci_max$max_al > thr),
                  round(thr)))
    }

    # ==========================================================
    # MODO A: boxplot e dots em faixas separadas
    # ==========================================================
    if (isTRUE(st$box_dots_separate)) {
      pdat[, locus_num := as.numeric(locus_label)]

      p <- ggplot(pdat, aes(y = allele2_est, fill = group, colour = group))

      if (isTRUE(st$show_violin)) {
        p <- p + geom_flat_violin(
          aes(x = locus_num + st$violin_nudge),
          trim   = FALSE,
          alpha  = 0.5,
          colour = NA,
          width  = st$violin_width
        )
      }

      p <- p + geom_point(
        aes(x = locus_num + st$dots_nudge),
        size        = st$point_size,
        alpha       = 0.6,
        position    = position_jitter(width = st$point_jitter, seed = 321),
        show.legend = FALSE
      )

      p <- p + geom_boxplot(
        aes(x = locus_num + st$box_nudge,
            group = interaction(locus_label, group)),
        width         = st$box_width,
        outlier.shape = NA,
        alpha         = 1,
        colour        = "black",
        position      = position_dodge(width = st$box_width),
        show.legend   = FALSE
      )

      p <- p +
        scale_fill_manual(values = group_colors, name = "Group") +
        scale_colour_manual(values = group_colors, guide = "none") +
        scale_x_continuous(
          breaks = seq_along(levels(pdat$locus_label)),
          labels = levels(pdat$locus_label)
        ) +
        scale_y_continuous(
          breaks = scales::breaks_pretty(n = 8),
          expand = expansion(mult = c(0.02, 0.05))
        ) +
        coord_flip() +
        labs(
          title    = title,
          subtitle = NULL,
          x        = NULL,
          y        = "Major-allele length (repeat units)",
          caption  = split_caption
        ) +
        theme_classic(base_size = 11) +
        theme(
          axis.ticks.y     = element_blank(),
          axis.line.y      = element_blank(),
          axis.title.y     = element_blank(),
          axis.text.y      = element_text(size = 10, color = "grey20"),
          axis.text.x      = element_text(size = 12, color = "grey20"),
          axis.title.x     = element_text(size = 11, face = "bold"),
          plot.title       = element_text(size = 13, hjust = 0.5, face = "bold"),
          plot.caption     = element_text(size = 10, color = "grey30",
                                          hjust = 0.5, face = "italic",
                                          margin = margin(t = 10)),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text       = element_text(size = 11, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_text(size = 10),
          panel.spacing    = unit(0.8, "lines"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)
        )

      if (use_facet) {
        p <- p + facet_wrap(~ size_group, scales = "free", ncol = 1)
      }

    # ==========================================================
    # MODO B: layout original
    # ==========================================================
    } else {
      p <- ggplot(pdat,
                  aes(x = locus_label, y = allele2_est,
                      fill = group, colour = group)) +
        geom_flat_violin(
          position = position_nudge(x = st$violin_nudge),
          trim     = FALSE,
          alpha    = 0.5,
          colour   = NA,
          width    = st$violin_width
        ) +
        geom_boxplot(
          width         = st$box_width,
          outlier.shape = NA,
          alpha         = 0.6,
          colour        = "black",
          position      = position_nudge(x = st$box_nudge),
          show.legend   = FALSE
        ) +
        geom_point(
          size        = st$point_size,
          alpha       = 0.5,
          position    = position_jitter(width = st$point_jitter, seed = 321),
          show.legend = FALSE
        ) +
        scale_fill_manual(values = group_colors, name = "Group") +
        scale_colour_manual(values = group_colors, guide = "none") +
        coord_flip() +
        scale_y_continuous(
          breaks = scales::breaks_pretty(n = 8),
          expand = expansion(mult = c(0.02, 0.05))
        ) +
        labs(
          title    = title,
          subtitle = NULL,
          x        = NULL,
          y        = "Major-allele length (repeat units)"
        ) +
        theme_classic(base_size = 11) +
        theme(
          axis.ticks.y     = element_blank(),
          axis.line.y      = element_blank(),
          axis.title.y     = element_blank(),
          axis.text.y      = element_text(size = 10, color = "grey20"),
          axis.text.x      = element_text(size = 9, color = "grey20"),
          axis.title.x     = element_text(size = 11, face = "bold"),
          plot.title       = element_text(size = 13, hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text       = element_text(size = 11, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_text(size = 10),
          panel.spacing    = unit(1.2, "lines"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)
        )
      # (facet_wrap(~ gse) removido — nao ha mais visao ALL)
    }

    out_file <- file.path(dest_dir,
      sprintf("%s_raincloud_per_locus_p%02d.png",
              gsub("[^A-Za-z0-9._-]", "_", intv), pg_i))

    n_rows_effective <- if (use_facet) n_loci * 1.6 else n_loci
    ggsave(
      filename  = out_file,
      plot      = p,
      width     = 12,
      height    = max(6, n_rows_effective * st$height_factor + st$height_offset),
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
# 2. Loop por intervenção — apenas por GSE
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

  # --- Apenas por GSE ---
  for (gse_i in sort(unique(dat$gse))) {
    dat_gse <- dat[gse == gse_i]
    cat(sprintf("  Gerando por GSE: %s\n", gse_i))
    write_patients_csv(dat_gse, intv, gse_i)
    make_raincloud(
      dat = dat_gse, intv = intv, gse_tag = gse_i,
      title = sprintf("DBSCAN outliers per locus - %s (%s)",
                      gsub("_", " ", intv), gse_i),
      subtitle = NULL
    )
  }
}

cat("\nConcluido.\n")