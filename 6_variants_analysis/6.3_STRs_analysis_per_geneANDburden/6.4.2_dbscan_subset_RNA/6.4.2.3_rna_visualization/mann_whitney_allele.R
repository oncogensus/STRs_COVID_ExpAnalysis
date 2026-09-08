#!/usr/bin/env Rscript
# mann_whitney_allele.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Teste de Mann-Whitney (Wilcoxon rank-sum) para comparar o tamanho dos
#   alelos (media de alelos e allele2_est) entre grupos case e control,
#   para cada STR localizado em genes DEGs do RNA-Seq.
#
#   Para CADA metrica (mean_allele e allele2), gera:
#     1) Tabela com estatisticas por STR (U, p, effect size, medias/grupos)
#     2) Manhattan plot (-log10 p) por STR
#     3) Boxplot para STRs significativos
#
#   Ao final, gera comparativo entre alelos (concordancia/discordancia).
#
# ENTRADAS
#   --str-catalog       STRs_analysis_dataset.tsv (per-sample x STR, com group)
#   --rna-gene-strs     rna_gene_strs.tsv (STRs que sao DEGs)
#   --out-dir           Diretorio de saida
#
# DEPENDENCIAS
#   R: data.table, ggplot2, scales, patchwork
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
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

path_str_catalog   <- parse_arg("--str-catalog")
path_rna_gene_strs <- parse_arg("--rna-gene-strs")
out_dir            <- parse_arg("--out-dir", ".")

missing <- c(
  "--str-catalog" = is.null(path_str_catalog),
  "--rna-gene-strs" = is.null(path_rna_gene_strs)
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
cat(sprintf("  %d linhas (sample x STR)\n", nrow(str_cat)))

cat("Carregando rna_gene_strs.tsv...\n")
rna_gene_strs <- fread(path_rna_gene_strs, header = TRUE, sep = "\t")
cat(sprintf("  %d pares gene x STR\n", nrow(rna_gene_strs)))

# ==========================================
# 2. Filter: only DEG STRs with valid group
# ==========================================
cat("\nFiltrando dados...\n")

de_strs <- unique(rna_gene_strs$strs_id)
str_deg <- str_cat[STRs_ID %in% de_strs]

str_deg <- str_deg[!is.na(group) & group != ""]
str_deg <- str_deg[group %in% c("case", "control")]

# Compute mean allele: (allele1_est + allele2_est) / 2
str_deg[, mean_allele := (allele1_est + allele2_est) / 2]

cat(sprintf("  STRs DEGs: %d unicos\n", length(unique(str_deg$STRs_ID))))
cat(sprintf("  Amostras: %d (case=%d, control=%d)\n",
            nrow(str_deg),
            sum(str_deg$group == "case"),
            sum(str_deg$group == "control")))

# ==========================================
# 3. Function: Mann-Whitney for one allele
# ==========================================
run_mann_whitney <- function(dat, allele_col, allele_label) {
  cat(sprintf("\n--- %s ---\n", allele_label))

  # Filter valid values
  d <- dat[!is.na(get(allele_col)) & get(allele_col) > 0]
  str_list <- unique(d$STRs_ID)
  cat(sprintf("  STRs com dado valido: %d\n", length(str_list)))

  results <- vector("list", length(str_list))
  for (i in seq_along(str_list)) {
    sid <- str_list[i]
    ds <- d[STRs_ID == sid]
    case_vals <- ds[group == "case", get(allele_col)]
    control_vals <- ds[group == "control", get(allele_col)]

    if (length(case_vals) < 2 || length(control_vals) < 2) {
      results[[i]] <- data.table(
        STRs_ID = sid, gene_name = ds$gene_name[1], region = ds$region[1],
        chrom = ds$chrom[1], start = ds$start[1],
        n_case = length(case_vals), n_control = length(control_vals),
        mean_case = mean(case_vals), mean_control = mean(control_vals),
        median_case = median(case_vals), median_control = median(control_vals),
        statistic_U = NA_real_, p_value = NA_real_, effect_size_r = NA_real_,
        mean_diff = mean(case_vals) - mean(control_vals),
        median_diff = median(case_vals) - median(control_vals),
        significant_bh = NA
      )
      next
    }

    wt <- wilcox.test(case_vals, control_vals, exact = FALSE, conf.int = TRUE)
    n1 <- length(case_vals)
    n2 <- length(control_vals)
    U <- as.numeric(wt$statistic)
    r <- 1 - (2 * U) / (n1 * n2)

    results[[i]] <- data.table(
      STRs_ID = sid, gene_name = ds$gene_name[1], region = ds$region[1],
      chrom = ds$chrom[1], start = ds$start[1],
      n_case = n1, n_control = n2,
      mean_case = mean(case_vals), mean_control = mean(control_vals),
      median_case = median(case_vals), median_control = median(control_vals),
      statistic_U = U, p_value = wt$p.value, effect_size_r = r,
      mean_diff = mean(case_vals) - mean(control_vals),
      median_diff = median(case_vals) - median(control_vals),
      significant_bh = NA
    )
  }

  res <- rbindlist(results)

  # BH correction
  valid_p <- !is.na(res$p_value)
  res[valid_p, p_adjusted := p.adjust(p_value, method = "BH")]
  res[, significant_bh := !is.na(p_adjusted) & p_adjusted < 0.05]
  res <- res[order(p_value)]

  n_tested <- sum(valid_p)
  n_sig <- sum(res$significant_bh, na.rm = TRUE)
  cat(sprintf("  Testes validos: %d\n", n_tested))
  cat(sprintf("  Significativos (BH<0.05): %d (%.1f%%)\n",
              n_sig, ifelse(n_tested > 0, n_sig / n_tested * 100, 0)))

  return(list(results = res, n_tested = n_tested, n_sig = n_sig))
}

# ==========================================
# 4. Run for both alleles
# ==========================================
res_mean <- run_mann_whitney(str_deg, "mean_allele", "Mean Allele ((A1+A2)/2)")
res_a2   <- run_mann_whitney(str_deg, "allele2_est", "Allele 2 (allele2_est)")

mw_mean <- res_mean$results
mw_a2   <- res_a2$results

# ==========================================
# 5. Save TSV for each allele
# ==========================================
out_mean <- file.path(out_dir, "mann_whitney_mean_allele_results.tsv")
out_a2   <- file.path(out_dir, "mann_whitney_allele2_results.tsv")
fwrite(mw_mean, out_mean, sep = "\t")
fwrite(mw_a2, out_a2, sep = "\t")
cat(sprintf("\nTabelas salvas:\n  %s (%d STRs)\n  %s (%d STRs)\n",
            out_mean, nrow(mw_mean), out_a2, nrow(mw_a2)))

# ==========================================
# 6. Manhattan plots (one per allele)
# ==========================================
cat("\nGerando Manhattan plots...\n")

make_manhattan <- function(mw, allele_label, n_tested, outfile) {
  pd <- mw[!is.na(p_value)]
  pd[, neg_log10_p := -log10(p_value)]
  pd <- pd[order(chrom, start)]
  pd[, x_pos := .I]
  pd[, chrom_num := as.integer(gsub("[^0-9]", "", chrom))]
  pd[is.na(chrom_num), chrom_num := 23]
  pd[, chrom_color := ifelse(chrom_num %% 2 == 0, "#4DAF4A", "#377EB8")]

  bonf_line <- -log10(0.05 / n_tested)
  nominal_line <- -log10(0.05)

  p <- ggplot(pd, aes(x = x_pos, y = neg_log10_p)) +
    geom_point(aes(color = chrom_color), size = 2, alpha = 0.8) +
    geom_hline(yintercept = bonf_line, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_hline(yintercept = nominal_line, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    scale_color_identity() +
    scale_x_continuous(
      breaks = pd[, .(x_pos = mean(x_pos)), by = chrom]$x_pos,
      labels = pd[, .(label = chrom[1]), by = chrom]$label,
      expand = c(0.02, 0)
    ) +
    labs(
      title = paste("Mann-Whitney U test:", allele_label, "case vs control"),
      subtitle = sprintf("Red = Bonferroni (p < %.2e); dotted = p < 0.05", 0.05 / n_tested),
      x = "Chromosome", y = expression(-log[10](p))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(outfile, p, width = 14, height = 6, dpi = 300, bg = "white")
  cat(sprintf("  %s\n", outfile))
}

make_manhattan(mw_mean, "Mean Allele", res_mean$n_tested,
               file.path(out_dir, "mann_whitney_manhattan_mean_allele.png"))
make_manhattan(mw_a2, "Allele 2", res_a2$n_tested,
               file.path(out_dir, "mann_whitney_manhattan_allele2.png"))

# ==========================================
# 7. Boxplots (one per allele)
# ==========================================
cat("\nGerando boxplots...\n")

make_boxplot <- function(mw, dat, allele_col, allele_label, outfile) {
  sig_strs <- mw[significant_bh == TRUE][order(p_value)]
  if (nrow(sig_strs) == 0) {
    cat(sprintf("  %s: nenhum STR significativo, pulando boxplot\n", allele_label))
    return(invisible(NULL))
  }
  if (nrow(sig_strs) > 12) sig_strs <- sig_strs[1:12]

  bd <- dat[STRs_ID %in% sig_strs$STRs_ID]
  bd <- bd[!is.na(get(allele_col)) & get(allele_col) > 0]
  sig_order <- sig_strs$STRs_ID
  bd[, STRs_ID := factor(STRs_ID, levels = sig_order)]

  p <- ggplot(bd, aes(x = STRs_ID, y = get(allele_col), fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(aes(color = group), width = 0.15, size = 0.8, alpha = 0.4) +
    scale_fill_manual(values = c("case" = "#E41A1C", "control" = "#377EB8")) +
    scale_color_manual(values = c("case" = "#E41A1C", "control" = "#377EB8")) +
    labs(
      title = paste(allele_label, ": case vs control (significant STRs)"),
      subtitle = "Mann-Whitney U test, BH-adjusted p < 0.05",
      x = "STR locus", y = paste(allele_label, "(repeat units)"),
      fill = "Group", color = "Group"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  ggsave(outfile, p, width = max(8, nrow(sig_strs) * 1.2), height = 6, dpi = 300, bg = "white")
  cat(sprintf("  %s\n", outfile))
}

make_boxplot(mw_mean, str_deg, "mean_allele", "Mean Allele",
             file.path(out_dir, "mann_whitney_boxplot_mean_allele.png"))
make_boxplot(mw_a2, str_deg, "allele2_est", "Allele 2",
             file.path(out_dir, "mann_whitney_boxplot_allele2.png"))

# ==========================================
# 8. Concordance between alleles
# ==========================================
cat("\nGerando comparativo entre alelos...\n")

# Merge results from both alleles
merged <- merge(
  mw_mean[, .(STRs_ID, gene_name, p_adjusted_mean = p_adjusted, sig_mean = significant_bh,
               effect_r_mean = effect_size_r, median_diff_mean = median_diff)],
  mw_a2[, .(STRs_ID, p_adjusted_a2 = p_adjusted, sig_a2 = significant_bh,
             effect_r_a2 = effect_size_r, median_diff_a2 = median_diff)],
  by = "STRs_ID", all = TRUE
)

merged[, concordance := fifelse(
  sig_mean == TRUE & sig_a2 == TRUE, "both",
  fifelse(sig_mean == TRUE, "mean_only",
  fifelse(sig_a2 == TRUE, "allele2_only", "none"))
)]

n_both <- sum(merged$concordance == "both", na.rm = TRUE)
n_mean_only <- sum(merged$concordance == "mean_only", na.rm = TRUE)
n_a2_only <- sum(merged$concordance == "allele2_only", na.rm = TRUE)
n_none <- sum(merged$concordance == "none", na.rm = TRUE)

cat(sprintf("  Concordancia entre metricas:\n"))
cat(sprintf("    Significativo em AMBOS:       %d\n", n_both))
cat(sprintf("    Significativo so mean_allele: %d\n", n_mean_only))
cat(sprintf("    Significativo so allele2:     %d\n", n_a2_only))
cat(sprintf("    Nao significativo em nenhum:  %d\n", n_none))

# Save concordance table
out_conc <- file.path(out_dir, "mann_whitney_concordance.tsv")
fwrite(merged, out_conc, sep = "\t")
cat(sprintf("\nTabela concordancia: %s\n", out_conc))

# Concordance plot
conc_data <- merged[!is.na(concordance)]
  conc_data[, concordance := factor(concordance,
                                   levels = c("both", "mean_only", "allele2_only", "none"))]

p_conc <- ggplot(conc_data, aes(x = effect_r_mean, y = effect_r_a2, color = concordance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(
    values = c("both" = "#E41A1C", "mean_only" = "#377EB8",
               "allele2_only" = "#4DAF4A", "none" = "grey70"),
    name = "Significance"
  ) +
  labs(
    title = "Effect size concordance: mean allele vs allele2",
    subtitle = "Rank-biserial correlation (r) for mean_allele vs allele2_est",
    x = expression(r[mean~allele]),
    y = expression(r[allele2])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

out_conc_plot <- file.path(out_dir, "mann_whitney_concordance_plot.png")
ggsave(out_conc_plot, p_conc, width = 8, height = 7, dpi = 300, bg = "white")
cat(sprintf("Plot concordancia: %s\n", out_conc_plot))

# ==========================================
# 9. Summary
# ==========================================
cat("\n", strrep("=", 55), "\n")
cat(" RESUMO: Mann-Whitney U test (case vs control)\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  Casos (case):              %d amostras\n",
            sum(str_deg$group == "case")))
cat(sprintf("  Controles (control):       %d amostras\n",
            sum(str_deg$group == "control")))
cat(sprintf("\n  {'':35s} {'Mean Allele':>12s} {'Allele 2':>12s}\n"))
cat(sprintf("  {'STRs testados':35s} %12d %12d\n", res_mean$n_tested, res_a2$n_tested))
cat(sprintf("  {'Significativos (BH<0.05)':35s} %12d %12d\n", res_mean$n_sig, res_a2$n_sig))
cat(sprintf("  {'%%':35s} %11.1f%% %11.1f%%\n",
            ifelse(res_mean$n_tested > 0, res_mean$n_sig / res_mean$n_tested * 100, 0),
            ifelse(res_a2$n_tested > 0, res_a2$n_sig / res_a2$n_tested * 100, 0)))
cat(sprintf("\n  Concordancia: %d em ambos, %d so mean_allele, %d so allele2\n",
            n_both, n_mean_only, n_a2_only))
cat(strrep("=", 55), "\n")
cat("Concluido.\n")
