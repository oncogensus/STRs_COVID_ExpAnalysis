#!/usr/bin/env Rscript
## plot_pathways.R  (ETAPA DE FIGURAS — le somente results/ + inputs)
## Gera as figuras de cross-validation a partir das tabelas produzidas por
## pathway_crossvalidation.R. Pode rodar independentemente (apos inspecionar os dados).
suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(UpSetR)
})

ROOT <- normalizePath("../../..")
p1_file    <- "../dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file    <- "../data/outlier_genes.txt"
resid_file <- "../dbscan_subset_GWAS/data/suggestive_strs_residuals.tsv"

ORIGINAL_8 <- c("ROBO2","ANK3","CDH12","NKAIN2","SEMA6D","KCNH1","KCNQ5","ST6GALNAC3")
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

## ======================================================================
## FIG A: dotplot de vias enriquecidas por pipeline (nominal p<0.05; cor = classe)
## ======================================================================
safe_read <- function(f) if (file.exists(f)) read_csv(f) else NULL

PIPE_LAB <- c("P1" = "GWAS-COVID filtered", "P2" = "Global outliers")

dfA <- bind_rows(
  safe_read("results/KEGG_P1_enrich.csv")    %>% mutate(source = "KEGG",    pipeline = "P1"),
  safe_read("results/KEGG_P2_enrich.csv")    %>% mutate(source = "KEGG",    pipeline = "P2"),
  safe_read("results/Reactome_P1_enrich.csv") %>% mutate(source = "Reactome", pipeline = "P1"),
  safe_read("results/Reactome_P2_enrich.csv") %>% mutate(source = "Reactome", pipeline = "P2")
) %>% filter(!is.na(pvalue) & pvalue < 0.05)

if (!is.null(dfA) && nrow(dfA) > 0) {
  cls <- tryCatch(read_delim("results/pathway_convergence_nominal_p005.tsv", delim = "\t"), error = function(e) NULL)
  if (!is.null(cls)) dfA <- left_join(dfA, cls %>% select(source, ID, class), by = c("source", "ID"))
  dfA <- dfA %>% mutate(n_gene = lengths(strsplit(geneID, "/")),
                        class = ifelse(is.na(.data$class), "None", as.character(.data$class)),
                        lab = paste(source, Description, sep = ": "))
  pA <- ggplot(dfA, aes(x = -log10(pvalue), y = reorder(lab, -log10(pvalue)))) +
    geom_point(aes(size = n_gene, color = class), alpha = 0.85) +
    scale_color_manual(values = c("High-Confidence" = "#2C3E50",
                                  "Hypothesis-Driven" = "#3498DB",
                                  "Agnostic-Specific" = "#c0392b", "None" = "#999999")) +
    facet_wrap(~ source, scales = "free_y") +
    theme_minimal() + labs(title = "Vias enriquecidas por pipeline (nominal p<0.05)",
                           x = "-log10(p bruto)", y = "", size = "n genes")
  ggsave("results/figA_dotplot.png", pA, width = 12, height = 10, dpi = 300, bg = "white")
  cat("Fig A ok\n")
} else cat("Sem vias significativas para Fig A (veja *_top_pvalue.csv).\n")

## ======================================================================
## FIG B: cnetplot das vias High-Confidence + heatmap de membership
## ======================================================================
cls_master  <- tryCatch(read_delim("results/pathway_convergence.tsv", delim = "\t"), error = function(e) NULL)
cls_nominal <- tryCatch(read_delim("results/pathway_convergence_nominal_p005.tsv", delim = "\t"), error = function(e) NULL)
if (!is.null(cls_nominal) && nrow(cls_nominal) > 0) {
  hc_k <- cls_nominal$ID[cls_nominal$class == "High-Confidence" & cls_nominal$source == "KEGG"]
  hc_r <- cls_nominal$ID[cls_nominal$class == "High-Confidence" & cls_nominal$source == "Reactome"]
  plots_B <- list()
  if (length(hc_k) > 0 && file.exists("results/KEGG_P1.rds")) {
    kk <- readRDS("results/KEGG_P1.rds")
    plots_B <- c(plots_B, list(cnetplot(kk, showCategory = hc_k, colorEdge = TRUE) +
                                 labs(title = "High-Confidence (KEGG)")))
  }
  if (length(hc_r) > 0 && file.exists("results/Reactome_P1.rds")) {
    rr <- readRDS("results/Reactome_P1.rds")
    plots_B <- c(plots_B, list(cnetplot(rr, showCategory = hc_r, colorEdge = TRUE) +
                                 labs(title = "High-Confidence (Reactome)")))
  }
  if (length(plots_B) > 0) {
    ggsave("results/figB_cnetplot.png", wrap_plots(plots_B, ncol = 1),
           width = 10, height = 6 * length(plots_B), dpi = 300, bg = "white")
    cat("Fig B cnetplot ok\n")
  }
}
if (!is.null(cls_master) && nrow(cls_master) > 0) {
  mem <- cls_master %>% filter(class != "None") %>%
    select(source, ID, Description, class, in_p1, in_p2) %>%
    pivot_longer(cols = c(in_p1, in_p2), names_to = "pip", values_to = "member") %>%
    mutate(pip = ifelse(pip == "in_p1", "P1", "P2"),
           lab = fct_reorder(paste(source, Description, sep = ": "), as.numeric(class)))
  pBmem <- ggplot(mem, aes(x = pip, y = lab, fill = member)) +
    geom_tile(color = "white") + scale_fill_manual(values = c("white", "#2C3E50")) +
    theme_minimal() + theme(axis.text.y = element_text(size = 7)) +
    labs(title = "Membership de vias por pipeline", x = "", y = "", fill = "na via")
  ggsave("results/figB_membership.png", pBmem, width = 6, height = 10, dpi = 300, bg = "white")
  cat("Fig B membership ok\n")
}

## ======================================================================
## FIG C: sobreposicao com hsa05171 (COVID-19) via msigdbr (offline)
## ======================================================================
p1_sym <- unique(toupper(str_trim(read_delim(p1_file, delim = "\t")$gene)))
p2_sym <- unique(toupper(str_trim(read_delim(p2_file, delim = "\t")$gene_name)))
tryCatch({
  library(msigdbr)
  ## msigdbr: curated pathways ficam em category = "C2", subcategory = "CP:KEGG"
  kegg_covid <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                                 subcategory = "CP:KEGG") %>% filter(gs_id == "hsa05171")
  covid_genes <- unique(kegg_covid$gene_symbol)
  dfc <- data.frame(gene = covid_genes, in_P1 = covid_genes %in% p1_sym,
                    in_P2 = covid_genes %in% p2_sym, stringsAsFactors = FALSE)
  write.table(dfc, "results/hsa05171_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  dfc_long <- dfc %>% pivot_longer(cols = c(in_P1, in_P2), names_to = "pip", values_to = "hit") %>%
    mutate(pip = ifelse(pip == "in_P1", "P1", "P2"), gene = fct_reorder(gene, hit))
  pC <- ggplot(dfc_long, aes(x = pip, y = gene, fill = hit)) + geom_tile(color = "white") +
    scale_fill_manual(values = c("white", "#c0392b")) + theme_minimal() +
    theme(axis.text.y = element_text(size = 7)) +
    labs(title = "Genes da via KEGG hsa05171 (COVID-19) presentes nos hits",
         x = "", y = "", fill = "no hit")
  ggsave("results/figC_hsa05171.png", pC, width = 5, height = 12, dpi = 300, bg = "white")
  cat("Fig C ok: hsa05171", sum(dfc$in_P1), "em P1,", sum(dfc$in_P2), "em P2\n")
}, error = function(e) cat("msigdbr indisponivel; pulando Fig C\n"))

## ======================================================================
## FIGURAS SUPLEMENTARES (gene-level DBSCAN replication)
## ======================================================================
if (file.exists(p1_file) && file.exists(resid_file)) {
  p1df <- read_delim(p1_file, delim = "\t")
  cv7 <- tryCatch({
    gc <- read_delim("results/gene_convergence.tsv", delim = "\t")
    gc$gene[gc$strategy == "Cross-validated"]
  }, error = function(e) character(0))
  p1df <- p1df %>% mutate(
    chrom_n = case_when(chrom == "chrX" ~ 23, chrom == "chrY" ~ 24,
                        TRUE ~ suppressWarnings(as.integer(str_replace(chrom, "chr", "")))),
    category = case_when(gene %in% cv7 ~ "Crossvalidated (7)",
                         gene %in% ORIGINAL_8 ~ "Original-8", TRUE ~ "Other"))
  p1df <- p1df %>% group_by(chrom_n) %>% mutate(x = gene_start / 1e6) %>% ungroup()
  p1df$x <- p1df$x + (as.numeric(factor(p1df$chrom_n)) - 1) * 300
  pMan <- ggplot(p1df, aes(x = x, y = -log10(as.numeric(best_p)), color = category)) +
    geom_point(size = 3) +
    geom_text(aes(label = ifelse(category != "Other", as.character(gene), "")), hjust = -0.1, size = 3) +
    scale_color_manual(values = c("Crossvalidated (7)" = "#c0392b",
                                  "Original-8" = "#f1c40f", "Other" = "#3498DB")) +
    theme_minimal() + labs(title = "Manhattan-like dos 20 genes com STR-outlier (DBSCAN)",
                           x = "Posicao cromossomica (Mb)", y = "-log10(p)")
  ggsave("results/fig_genelevel_manhattan.png", pMan, width = 10, height = 6, dpi = 300, bg = "white")

  resid <- read_delim(resid_file, delim = "\t")
  out7 <- p1df %>% filter(gene %in% cv7) %>% select(gene, strs_id, outlier_samples, outlier_residuals)
  plot_data <- out7 %>% left_join(resid, by = c("strs_id" = "STRs_ID")) %>%
    mutate(is_out = sample_id == outlier_samples)
  pRes <- ggplot(plot_data, aes(x = sample_id, y = allele2_residuals)) +
    geom_point(aes(color = is_out), size = 1.2) + facet_wrap(~ gene, scales = "free_x") +
    geom_hline(aes(yintercept = as.numeric(outlier_residuals)), linetype = 2, color = "#c0392b") +
    scale_color_manual(values = c("FALSE" = "#999999", "TRUE" = "#c0392b")) +
    theme_minimal() + theme(axis.text.x = element_blank()) +
    labs(title = "Residuos DBSCAN (subset COVID) - 7 genes crossvalidados",
         x = "", y = "allele2_residuals", color = "outlier")
  ggsave("results/fig_genelevel_residual.png", pRes, width = 12, height = 8, dpi = 300, bg = "white")
  cat("Figuras gene-level ok\n")
} else cat("Inputs gene-level ausentes; pulando suplementares.\n")

## ======================================================================
## FIG D: UpSet — sobreposicao P1 x P2 (genes com STR-outlier)
## ======================================================================
tryCatch({
  upset_lst <- list(`GWAS-filtered` = p1_sym, `Agnostic` = p2_sym)
  png("results/fig_gene_upset.png", width = 8, height = 6, units = "in", res = 300, bg = "white")
  upset(fromList(upset_lst), sets = c("GWAS-filtered", "Agnostic"),
        order.by = "freq", number.angles = 30,
        text.scale = c(1.4, 1.2, 1.2, 1, 1.4, 1.2),
        mainbar.y.label = "Genes por intersecao", sets.x.label = "Total por estrategia")
  dev.off()
  cat("Fig D UpSet ok\n")
}, error = function(e) cat("UpSet falhou:", conditionMessage(e), "\n"))

cat("\n=== Pronto: figuras em results/ ===\n")
