#!/usr/bin/env Rscript
## pathway_crossvalidation.R
## Cross-validation funcional de vias (KEGG + Reactome) entre dois pipelines de STR-outliers:
##   P1 (guiado / Tier 1+3): outliers do DBSCAN restritos a genes GWAS-suggestivos COVID (p<1e-5)
##   P2 (agnostico / genome-wide): outliers do DBSCAN detectados em todo o genoma
## Classificacao das vias emergentes:
##   High-Confidence / Robust Signal : via significativa (FDR<0.05) em P1 E P2
##   Hypothesis-Driven Target        : via significativa so em P1 (ruido de FDR reduzido pelo filtro GWAS)
##   Agnostic-Specific               : via significativa so em P2
##
## Adaptado do padrao de KEEG_analysis/KEEG_df_import.r (mapeamento SYMBOL->ENTREZ, enrichKEGG,
## setReadable, saidas .csv/.rds). Diferencas: (1) nossos genes sao SYMBOL (nao ENSEMBL);
## (2) dois pipelines; (3) acrescenta ReactomePA e as figuras de cross-validation.
##
## Uso: Rscript pathway_crossvalidation.R   (roda no env micromamba r_enrich_env)
suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(readr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(ggplot2)
  library(patchwork)
})

ROOT <- normalizePath("../../..")   # STRs_COVID_Analysis
p1_file     <- "../dbscan_subset/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file     <- "../litcovid_validation/data/outlier_genes.txt"
catalog     <- file.path(ROOT, "samples/STRs_analysis_dataset.tsv")
resid_file  <- "../dbscan_subset/data/suggestive_strs_residuals.tsv"
cross_file  <- "../dbscan_subset/results/crossvalidated_genes.tsv"

ORIGINAL_8 <- c("ROBO2","ANK3","CDH12","NKAIN2","SEMA6D","KCNH1","KCNQ5","ST6GALNAC3")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

## ======================================================================
## 1. CARREGAR LISTAS DE GENES (SIMBOLOS)
## ======================================================================
p1 <- read_delim(p1_file, delim = "\t")
p1_sym <- unique(toupper(str_trim(p1$gene)))
cat("[P1] genes sugestivos com STR-outlier:", length(p1_sym), "\n")

p2 <- read_delim(p2_file, delim = "\t")
p2_sym <- unique(toupper(str_trim(p2$gene_name)))
cat("[P2] genes genome-wide com STR-outlier:", length(p2_sym), "\n")

cat_tbl <- read_delim(catalog, delim = "\t")
gene_col <- names(cat_tbl)[grep("gene_name|\\bgene\\b", tolower(names(cat_tbl)))[1]]
bg_sym <- unique(toupper(str_trim(cat_tbl[[gene_col]])))
bg_sym <- bg_sym[!is.na(bg_sym) & bg_sym != ""]
cat("[Background] genes com STR na coorte:", length(bg_sym), "\n")

## ======================================================================
## 2. MAPEAMENTO SYMBOL -> ENTREZ (org.Hs.eg.db)
## ======================================================================
sym2entrez <- function(syms) {
  mapIds(org.Hs.eg.db, keys = syms, column = "ENTREZID",
         keytype = "SYMBOL", multiVals = "first")
}
bg_entrez <- sym2entrez(bg_sym); bg_entrez <- bg_entrez[!is.na(bg_entrez)]
cat("[Background] mapeados para ENTREZ:", length(bg_entrez), "\n")

p1_mapped <- sym2entrez(p1_sym); p1_mapped <- p1_mapped[!is.na(p1_mapped)]
p2_mapped <- sym2entrez(p2_sym); p2_mapped <- p2_mapped[!is.na(p2_mapped)]
p1_entrez <- intersect(p1_mapped, bg_entrez)
p2_entrez <- intersect(p2_mapped, bg_entrez)
cat("[P1] foreground ENTREZ:", length(p1_entrez),
    " | [P2] foreground ENTREZ:", length(p2_entrez), "\n")

if (length(p1_entrez) < 5) stop("P1 muito pequeno para enrichment")
if (length(p2_entrez) < 5) stop("P2 muito pequeno para enrichment")

## ======================================================================
## 3. ENRICHMENT (KEGG + Reactome), FDR<0.05 (BH)
## ======================================================================
run_kegg <- function(genes, bg) {
  res <- enrichKEGG(gene = genes, organism = "hsa", keyType = "ncbi-geneid",
                    universe = bg, pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", qvalueCutoff = 1.0)
  if (!is.null(res) && nrow(res@result) > 0)
    res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  res
}
run_reactome <- function(genes, bg) {
  res <- tryCatch({
    if (!requireNamespace("msigdbr", quietly = TRUE))
      stop("msigdbr indisponivel")
    rs <- msigdbr(species = "Homo sapiens", category = "CP", subcategory = "REACTOME") %>%
      dplyr::select(gs_id, gs_name, entrez_gene)
    rs <- rs[!is.na(rs$entrez_gene) & rs$entrez_gene != "", ]
    t2g <- unique(rs %>% dplyr::select(gs_id, entrez_gene))
    t2n <- unique(rs %>% dplyr::select(gs_id, gs_name))
    er <- enricher(gene = genes, universe = bg, pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", qvalueCutoff = 1.0,
                   TERM2GENE = t2g, TERM2NAME = t2n)
    if (!is.null(er) && nrow(er@result) > 0)
      er <- setReadable(er, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    er
  }, error = function(e) { cat("Reactome (msigdbr) falhou:", conditionMessage(e), "\n"); NULL })
  res
}

kk_p1 <- run_kegg(p1_entrez, bg_entrez)
kk_p2 <- run_kegg(p2_entrez, bg_entrez)
re_p1 <- run_reactome(p1_entrez, bg_entrez)
re_p2 <- run_reactome(p2_entrez, bg_entrez)

extract_sig <- function(res, source_name, pipeline) {
  if (is.null(res) || nrow(res@result) == 0) return(data.frame())
  r <- as.data.frame(res@result)
  r <- r[!is.na(r$p.adjust) & r$p.adjust < 0.05, ]
  if (nrow(r) == 0) return(data.frame())
  data.frame(source = source_name, pipeline = pipeline,
             ID = r$ID, Description = r$Description,
             pvalue = r$pvalue, p.adjust = r$p.adjust,
             geneID = r$geneID, stringsAsFactors = FALSE)
}
sig <- bind_rows(
  extract_sig(kk_p1, "KEGG",    "P1"),
  extract_sig(kk_p2, "KEGG",    "P2"),
  extract_sig(re_p1, "Reactome", "P1"),
  extract_sig(re_p2, "Reactome", "P2")
)
cat("Vias significativas (FDR<0.05):", nrow(sig), "\n")

## salvar resultados individuais (padrao do base)
save_individual <- function(res, tag) {
  if (!is.null(res) && nrow(res@result) > 0) {
    saveRDS(res, sprintf("results/%s.rds", tag))
    write.csv(as.data.frame(res), sprintf("results/%s_enrich.csv", tag), row.names = FALSE)
  }
}
save_individual(kk_p1, "KEGG_P1"); save_individual(kk_p2, "KEGG_P2")
save_individual(re_p1, "Reactome_P1"); save_individual(re_p2, "Reactome_P2")

## ======================================================================
## 4. CLASSIFICACAO DE CROSS-VALIDATION (nivel via)
## ======================================================================
cls <- sig %>% group_by(source, ID) %>% summarise(
  Description = first(Description),
  p1 = min(p.adjust[pipeline == "P1"], na.rm = TRUE),
  p2 = min(p.adjust[pipeline == "P2"], na.rm = TRUE),
  in_p1 = any(pipeline == "P1"),
  in_p2 = any(pipeline == "P2"), .groups = "drop") %>% mutate(
  class = case_when(in_p1 & in_p2 ~ "High-Confidence",
                    in_p1            ~ "Hypothesis-Driven",
                    in_p2            ~ "Agnostic-Specific",
                    TRUE             ~ "None"))
cls$class <- factor(cls$class,
                    levels = c("High-Confidence","Hypothesis-Driven","Agnostic-Specific"))
write.table(cls, "results/pathway_convergence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("High-Confidence:", sum(cls$class == "High-Confidence"),
    "| Hypothesis-Driven:", sum(cls$class == "Hypothesis-Driven"),
    "| Agnostic-Specific:", sum(cls$class == "Agnostic-Specific"), "\n")

## ======================================================================
## 5. FIGURAS
## ======================================================================
## Fig A: dotplot de vias enriquecidas por pipeline (cor = classe)
sig <- sig %>% mutate(n_gene = lengths(strsplit(geneID, "/")))
sig <- left_join(sig, cls %>% select(source, ID, class), by = c("source", "ID"))
pA <- ggplot(sig, aes(x = -log10(p.adjust),
                      y = reorder(paste(source, Description, sep = ": "), -log10(p.adjust)))) +
  geom_point(aes(size = n_gene, color = class), alpha = 0.85) +
  scale_color_manual(values = c("High-Confidence" = "#2C3E50",
                                "Hypothesis-Driven" = "#3498DB",
                                "Agnostic-Specific" = "#c0392b")) +
  facet_wrap(~ source, scales = "free_y") +
  theme_minimal() + labs(title = "Vias enriquecidas por pipeline (FDR<0.05)",
                         x = "-log10(FDR)", y = "", size = "n genes")
ggsave("results/figA_dotplot.png", pA, width = 12, height = 10, dpi = 300, bg = "white")

## Fig B: cnetplot das vias High-Confidence (reusa padrao do base) + heatmap de membership
hc_k <- cls$ID[cls$class == "High-Confidence" & cls$source == "KEGG"]
hc_r <- cls$ID[cls$class == "High-Confidence" & cls$source == "Reactome"]
plots_B <- list()
if (length(hc_k) > 0 && !is.null(kk_p1)) {
  plots_B <- c(plots_B, list(cnetplot(kk_p1, showCategory = hc_k, colorEdge = TRUE) +
                                labs(title = "High-Confidence (KEGG)")))
}
if (length(hc_r) > 0 && !is.null(re_p1)) {
  plots_B <- c(plots_B, list(cnetplot(re_p1, showCategory = hc_r, colorEdge = TRUE) +
                                labs(title = "High-Confidence (Reactome)")))
}
if (length(plots_B) > 0) {
  (wrap_plots(plots_B, ncol = 1)) +
    ggsave("results/figB_cnetplot.png", width = 10, height = 6 * length(plots_B), dpi = 300, bg = "white")
}
## heatmap de membership (sem pacotes extras)
mem <- cls %>% filter(class != "None") %>%
  mutate(pipeline = ifelse(in_p1, "P1", NA), pipeline2 = ifelse(in_p2, "P2", NA)) %>%
  select(source, ID, Description, class, in_p1, in_p2) %>%
  pivot_longer(cols = c(in_p1, in_p2), names_to = "pip", values_to = "member") %>%
  mutate(pip = ifelse(pip == "in_p1", "P1", "P2"),
         lab = paste(source, Description, sep = ": ")) %>%
  mutate(lab = fct_reorder(lab, as.numeric(class)))
pBmem <- ggplot(mem, aes(x = pip, y = lab, fill = member)) +
  geom_tile(color = "white") + scale_fill_manual(values = c("white", "#2C3E50")) +
  theme_minimal() + theme(axis.text.y = element_text(size = 7)) +
  labs(title = "Membership de vias por pipeline", x = "", y = "", fill = "na via")
ggsave("results/figB_membership.png", pBmem, width = 6, height = 10, dpi = 300, bg = "white")

## Fig C: sobreposicao com hsa05171 (COVID-19) via msigdbr (offline)
tryCatch({
  library(msigdbr)
  kegg_covid <- msigdbr(species = "Homo sapiens", category = "CP", subcategory = "KEGG") %>%
    filter(gs_id == "hsa05171")
  covid_genes <- unique(kegg_covid$gene_symbol)
  dfc <- data.frame(gene = covid_genes,
                    in_P1 = covid_genes %in% p1_sym,
                    in_P2 = covid_genes %in% p2_sym, stringsAsFactors = FALSE)
  write.table(dfc, "results/hsa05171_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  dfc_long <- dfc %>% pivot_longer(cols = c(in_P1, in_P2), names_to = "pip", values_to = "hit") %>%
    mutate(pip = ifelse(pip == "in_P1", "P1", "P2"),
           gene = fct_reorder(gene, hit))
  pC <- ggplot(dfc_long, aes(x = pip, y = gene, fill = hit)) + geom_tile(color = "white") +
    scale_fill_manual(values = c("white", "#c0392b")) + theme_minimal() +
    theme(axis.text.y = element_text(size = 7)) +
    labs(title = "Genes da via KEGG hsa05171 (COVID-19) presentes nos hits",
         x = "", y = "", fill = "no hit")
  ggsave("results/figC_hsa05171.png", pC, width = 5, height = 12, dpi = 300, bg = "white")
  cat("hsa05171: ", sum(dfc$in_P1), " genes em P1, ", sum(dfc$in_P2), " em P2\n")
}, error = function(e) cat("msigdbr indisponivel; pulando Fig C\n"))

## ======================================================================
## 6. FIGURAS SUPLEMENTARES (gene-level DBSCAN replication)
## ======================================================================
p1df <- read_delim(p1_file, delim = "\t")
cv7 <- unique(read_delim(cross_file, delim = "\t")$gene)
p1df <- p1df %>% mutate(
  chrom_n = case_when(chrom == "chrX" ~ 23, chrom == "chrY" ~ 24,
                      TRUE ~ suppressWarnings(as.integer(str_replace(chrom, "chr", "")))),
  category = case_when(gene %in% cv7 ~ "Crossvalidated (7)",
                       gene %in% ORIGINAL_8 ~ "Original-8",
                       TRUE ~ "Other"))
## posicao cumulativa grosseira (20 pontos)
p1df <- p1df %>% group_by(chrom_n) %>% mutate(x = gene_start / 1e6) %>% ungroup()
p1df$x <- p1df$x + (as.numeric(factor(p1df$chrom_n)) - 1) * 300
pMan <- ggplot(p1df, aes(x = x, y = -log10(as.numeric(best_p)), color = category)) +
  geom_point(size = 3) +
  geom_text(aes(label = ifelse(category != "Other", as.character(gene), "")),
            hjust = -0.1, size = 3) +
  scale_color_manual(values = c("Crossvalidated (7)" = "#c0392b",
                                "Original-8" = "#f1c40f", "Other" = "#3498DB")) +
  theme_minimal() + labs(title = "Manhattan-like dos 20 genes com STR-outlier (DBSCAN)",
                         x = "Posicao cromossomica (Mb)", y = "-log10(p)")
ggsave("results/fig_genelevel_manhattan.png", pMan, width = 10, height = 6, dpi = 300, bg = "white")

## scatter de residuos DBSCAN facetado nos 7 crossvalidados
resid <- read_delim(resid_file, delim = "\t")   # STRs_ID, allele2_residuals, sample_id, group
out7 <- p1df %>% filter(gene %in% cv7) %>%
  select(gene, strs_id, outlier_samples, outlier_residuals)
plot_data <- out7 %>% left_join(resid, by = c("strs_id" = "STRs_ID"))
plot_data <- plot_data %>% mutate(is_out = sample_id == outlier_samples)
pRes <- ggplot(plot_data, aes(x = sample_id, y = allele2_residuals)) +
  geom_point(aes(color = is_out), size = 1.2) +
  facet_wrap(~ gene, scales = "free_x") +
  geom_hline(aes(yintercept = as.numeric(outlier_residuals)), linetype = 2, color = "#c0392b") +
  scale_color_manual(values = c("FALSE" = "#999999", "TRUE" = "#c0392b")) +
  theme_minimal() + theme(axis.text.x = element_blank()) +
  labs(title = "Residuos DBSCAN (subset COVID) - 7 genes crossvalidados",
       x = "", y = "allele2_residuals", color = "outlier")
ggsave("results/fig_genelevel_residual.png", pRes, width = 12, height = 8, dpi = 300, bg = "white")

cat("\n=== Pronto: resultados em results/ ===\n")
