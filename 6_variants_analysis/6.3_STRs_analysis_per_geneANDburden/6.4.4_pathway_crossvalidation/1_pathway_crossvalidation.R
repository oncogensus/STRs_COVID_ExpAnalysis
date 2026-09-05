#!/usr/bin/env Rscript
## pathway_crossvalidation.R  (ETAPA DE DADOS — so tabelas)
## Cross-validation funcional de vias (KEGG + Reactome) entre dois pipelines de STR-outliers:
##   P1 (guiado): outliers DBSCAN restritos a genes GWAS-suggestivos COVID (p<1e-5)
##   P2 (agnostico): outliers DBSCAN genome-wide
## Gera TABELAS em results/ (enrichment + classificacao). As FIGURAS estao em plot_pathways.R.
suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(AnnotationDbi)
})

ROOT <- normalizePath("../../../..")
p1_file <- "../../6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file <- "../data/outlier_genes.txt"
catalog <- file.path(ROOT, "samples/STRs_analysis_dataset.tsv")

## ======================================================================
## 0. ARGUMENTOS DE LINHA DE COMANDO (parametrizacao GWAS x RNA)
##     Rscript 1_pathway_crossvalidation.R [--p1-file P] [--p2-file P]
##                                         [--p1-label L] [--p2-label L]
##                                         [--out BASE]
##   Exemplo (comparacao GWAS vs RNA):
##     Rscript 1_pathway_crossvalidation.R \
##       --p1-file ../../6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv \
##       --p2-file ../../6.4.2_dbscan_subset_RNA/6.4.2.2_RNA_matrix/results/rna_outlier_genes.tsv \
##       --p1-label GWAS --p2-label RNA \
##       --out pathway_convergence_gwas_rna
##   Sem argumentos: comportamento retroativo (P1=GWAS, P2=agnostico).
## ======================================================================
get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
args <- commandArgs(trailingOnly = TRUE)
p1_file <- get_opt(args, "--p1-file", p1_file)
p2_file <- get_opt(args, "--p2-file", p2_file)
p1_lab  <- get_opt(args, "--p1-label", "P1")
p2_lab  <- get_opt(args, "--p2-label", "P2")
out_tag <- get_opt(args, "--out", "pathway_convergence")

cat("[args] p1_file =", p1_file, "\n")
cat("[args] p2_file =", p2_file, "\n")
cat("[args] labels  =", p1_lab, "vs", p2_lab, "\n")
cat("[args] out tag =", out_tag, "\n")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

## ======================================================================
## 1. LISTAS DE GENES (SIMBOLOS)
## ======================================================================
read_gene_list <- function(path) {
  d <- read_delim(path, delim = "\t")
  gc <- grep("^gene(_name)?$", colnames(d), value = TRUE)[1]
  if (is.na(gc)) stop("coluna 'gene'/'gene_name' nao encontrada em: ", path)
  unique(toupper(str_trim(d[[gc]])))
}

p1_sym <- read_gene_list(p1_file)
cat(sprintf("[%s] genes sugestivos com STR-outlier: %d\n", p1_lab, length(p1_sym)))

p2_sym <- read_gene_list(p2_file)
cat(sprintf("[%s] genes com STR-outlier: %d\n", p2_lab, length(p2_sym)))

cat_tbl <- read_delim(catalog, delim = "\t")
gene_col <- names(cat_tbl)[grep("gene_name|\\bgene\\b", tolower(names(cat_tbl)))[1]]
bg_sym <- unique(toupper(str_trim(cat_tbl[[gene_col]])))
bg_sym <- bg_sym[!is.na(bg_sym) & bg_sym != ""]
cat("[Background] genes com STR na coorte:", length(bg_sym), "\n")

## ======================================================================
## 2. MAPEAMENTO SYMBOL -> ENTREZ
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
cat(sprintf("[%s] foreground ENTREZ: %d | [%s] foreground ENTREZ: %d\n",
            p1_lab, length(p1_entrez), p2_lab, length(p2_entrez)))

if (length(p1_entrez) < 5) stop("P1 muito pequeno para enrichment")
if (length(p2_entrez) < 5) stop("P2 muito pequeno para enrichment")

## ======================================================================
## 3. ENRICHMENT (KEGG online + Reactome via msigdbr/enricher, FDR<0.05)
## ======================================================================
## Retorna TODAS as vias testadas (cutoff=1) para re-filtrar por threshold em R.
run_kegg <- function(genes, bg) {
  res <- enrichKEGG(gene = genes, organism = "hsa", keyType = "ncbi-geneid",
                    universe = bg, pvalueCutoff = 1.0,
                    pAdjustMethod = "BH", qvalueCutoff = 1.0)
  if (!is.null(res) && nrow(res@result) > 0)
    res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  res
}
run_reactome <- function(genes, bg) {
  res <- tryCatch({
    if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr indisponivel")
    ## msigdbr: curated pathways ficam em category = "C2", subcategory = "CP:REACTOME"
    rs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                           subcategory = "CP:REACTOME") %>%
      dplyr::select(gs_id, gs_name, entrez_gene)
    rs <- rs[!is.na(rs$entrez_gene) & rs$entrez_gene != "", ]
    t2g <- unique(rs %>% dplyr::select(gs_id, entrez_gene))
    t2n <- unique(rs %>% dplyr::select(gs_id, gs_name))
    er <- enricher(gene = genes, universe = bg, pvalueCutoff = 1.0,
                   pAdjustMethod = "BH", qvalueCutoff = 1.0,
                   TERM2GENE = t2g, TERM2NAME = t2n)
    if (!is.null(er) && nrow(er@result) > 0)
      er <- setReadable(er, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    er
  }, error = function(e) { cat("Reactome (msigdbr) falhou:", conditionMessage(e), "\n"); NULL })
}

kk_p1 <- run_kegg(p1_entrez, bg_entrez)
kk_p2 <- run_kegg(p2_entrez, bg_entrez)
re_p1 <- run_reactome(p1_entrez, bg_entrez)
re_p2 <- run_reactome(p2_entrez, bg_entrez)

## Extrai TODAS as vias testadas (com p bruto e FDR) para re-filtrar por threshold.
extract_all <- function(res, source_name, pipeline) {
  if (is.null(res) || nrow(res@result) == 0) return(data.frame())
  r <- as.data.frame(res@result)
  data.frame(source = source_name, pipeline = pipeline,
             ID = r$ID, Description = r$Description,
             pvalue = r$pvalue, p.adjust = r$p.adjust,
             geneID = r$geneID, stringsAsFactors = FALSE)
}
all_sig <- bind_rows(
  extract_all(kk_p1, "KEGG",     p1_lab),
  extract_all(kk_p2, "KEGG",     p2_lab),
  extract_all(re_p1, "Reactome", p1_lab),
  extract_all(re_p2, "Reactome", p2_lab))
cat("[all_sig] vias testadas retidas:", nrow(all_sig), "\n")

## debug: top pathways por pvalue (mesmo que nenhuma passe FDR) para inspecao
dump_top <- function(res, tag) {
  if (!is.null(res) && nrow(res@result) > 0) {
    d <- as.data.frame(res@result) %>% arrange(pvalue)
    write.csv(head(d, 15), sprintf("results/%s_top_pvalue.csv", tag), row.names = FALSE)
  }
}
dump_top(kk_p1, paste0("KEGG_", p1_lab)); dump_top(kk_p2, paste0("KEGG_", p2_lab))
dump_top(re_p1, paste0("Reactome_", p1_lab)); dump_top(re_p2, paste0("Reactome_", p2_lab))

save_individual <- function(res, tag) {
  if (!is.null(res) && nrow(res@result) > 0) {
    saveRDS(res, sprintf("results/%s.rds", tag))
    write.csv(as.data.frame(res), sprintf("results/%s_enrich.csv", tag), row.names = FALSE)
  }
}
save_individual(kk_p1, paste0("KEGG_", p1_lab)); save_individual(kk_p2, paste0("KEGG_", p2_lab))
save_individual(re_p1, paste0("Reactome_", p1_lab)); save_individual(re_p2, paste0("Reactome_", p2_lab))

## ======================================================================
## 4. CLASSIFICACAO DE CROSS-VALIDATION (nivel via) POR THRESHOLD
##    - nominal: p bruto < 0.05
##    - fdr005 : FDR (BH)  < 0.05
##    - fdr001 : FDR (BH)  < 0.01
## Cada threshold gera seu proprio CSV/TSV (pedido do usuario).
## ======================================================================
minna <- function(x) if (length(x) == 0) NA else min(x, na.rm = TRUE)

summarise_cls <- function(sig, tag) {
  if (nrow(sig) == 0) {
    cat(sprintf("AVISO: nenhuma via em '%s'. Veja results/*_top_pvalue.csv.\n", tag))
    cls <- data.frame(source = character(), ID = character(), Description = character(),
                      p1_pvalue = numeric(), p2_pvalue = numeric(),
                      p1_p.adjust = numeric(), p2_p.adjust = numeric(),
                      in_p1 = logical(), in_p2 = logical(),
                      class = character(), stringsAsFactors = FALSE)
  } else {
    cls <- sig %>% group_by(source, ID) %>% summarise(
      Description  = Description[1],
      p1_pvalue    = minna(pvalue[pipeline == p1_lab]),
      p2_pvalue    = minna(pvalue[pipeline == p2_lab]),
      p1_p.adjust  = minna(p.adjust[pipeline == p1_lab]),
      p2_p.adjust  = minna(p.adjust[pipeline == p2_lab]),
      in_p1 = any(pipeline == p1_lab),
      in_p2 = any(pipeline == p2_lab), .groups = "drop") %>% mutate(
      class = case_when(in_p1 & in_p2 ~ "High-Confidence",
                        in_p1            ~ "Hypothesis-Driven",
                        in_p2            ~ "Agnostic-Specific",
                        TRUE             ~ "None"))
  }
  write.table(cls, sprintf("results/%s_%s.tsv", out_tag, tag),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(cls, sprintf("results/%s_%s.csv", out_tag, tag), row.names = FALSE)
  if (nrow(cls) > 0)
    cat(sprintf("[%s] High-Confidence: %d | Hypothesis-Driven: %d | Agnostic-Specific: %d\n",
                tag, sum(cls$class == "High-Confidence"),
                sum(cls$class == "Hypothesis-Driven"),
                sum(cls$class == "Agnostic-Specific")))
  invisible(cls)
}

sig_nominal <- all_sig %>% filter(!is.na(pvalue) & pvalue < 0.05)
sig_fdr005 <- all_sig %>% filter(!is.na(p.adjust) & p.adjust < 0.05)
sig_fdr001 <- all_sig %>% filter(!is.na(p.adjust) & p.adjust < 0.01)

cat("=== Cross-validation por threshold ===\n")
build_convergence_nominal <- summarise_cls(sig_nominal, "nominal_p005")
build_convergence_fdr005  <- summarise_cls(sig_fdr005,  "fdr005")
build_convergence_fdr001  <- summarise_cls(sig_fdr001,  "fdr001")

## ======================================================================
## 5. TABELA-MAE unindo os tres thresholds (p bruto + FDR 0.1 + FDR 0.5)
## ======================================================================
master <- all_sig %>% group_by(source, ID) %>% summarise(
  Description  = Description[1],
  p1_pvalue    = minna(pvalue[pipeline == p1_lab]),
  p2_pvalue    = minna(pvalue[pipeline == p2_lab]),
  p1_p.adjust  = minna(p.adjust[pipeline == p1_lab]),
  p2_p.adjust  = minna(p.adjust[pipeline == p2_lab]),
  in_p1 = any(pipeline == p1_lab),
  in_p2 = any(pipeline == p2_lab), .groups = "drop") %>% mutate(
  pass_nominal = (in_p1 & !is.na(p1_pvalue) & p1_pvalue < 0.05) |
                 (in_p2 & !is.na(p2_pvalue) & p2_pvalue < 0.05),
  pass_fdr005  = (in_p1 & !is.na(p1_p.adjust) & p1_p.adjust < 0.05) |
                 (in_p2 & !is.na(p2_p.adjust) & p2_p.adjust < 0.05),
  pass_fdr001  = (in_p1 & !is.na(p1_p.adjust) & p1_p.adjust < 0.01) |
                 (in_p2 & !is.na(p2_p.adjust) & p2_p.adjust < 0.01),
  class = case_when(in_p1 & in_p2 ~ "High-Confidence",
                    in_p1            ~ "Hypothesis-Driven",
                    in_p2            ~ "Agnostic-Specific",
                    TRUE             ~ "None"))
write.table(master, sprintf("results/%s.tsv", out_tag), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Tabela-mae: results/%s.tsv (%d vias unicas)\n", out_tag, nrow(master)))

cat("\n=== Pronto: tabelas em results/ (figuras em plot_pathways.R) ===\n")
