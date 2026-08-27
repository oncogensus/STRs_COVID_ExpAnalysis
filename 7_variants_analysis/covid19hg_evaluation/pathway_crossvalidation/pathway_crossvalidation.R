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

ROOT <- normalizePath("../../..")
p1_file <- "../dbscan_subset/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file <- "../litcovid_validation/data/outlier_genes.txt"
catalog <- file.path(ROOT, "samples/STRs_analysis_dataset.tsv")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

## ======================================================================
## 1. LISTAS DE GENES (SIMBOLOS)
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
cat("[P1] foreground ENTREZ:", length(p1_entrez),
    " | [P2] foreground ENTREZ:", length(p2_entrez), "\n")

if (length(p1_entrez) < 5) stop("P1 muito pequeno para enrichment")
if (length(p2_entrez) < 5) stop("P2 muito pequeno para enrichment")

## ======================================================================
## 3. ENRICHMENT (KEGG online + Reactome via msigdbr/enricher, FDR<0.05)
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
    if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr indisponivel")
    rs <- msigdbr::msigdbr(species = "Homo sapiens", category = "CP",
                           subcategory = "REACTOME") %>%
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
  extract_sig(re_p2, "Reactome", "P2"))
cat("Vias significativas (FDR<0.05):", nrow(sig), "\n")

## debug: top pathways por pvalue (mesmo que nenhuma passe FDR) para inspecao
dump_top <- function(res, tag) {
  if (!is.null(res) && nrow(res@result) > 0) {
    d <- as.data.frame(res@result) %>% arrange(pvalue)
    write.csv(head(d, 15), sprintf("results/%s_top_pvalue.csv", tag), row.names = FALSE)
  }
}
dump_top(kk_p1, "KEGG_P1"); dump_top(kk_p2, "KEGG_P2")
dump_top(re_p1, "Reactome_P1"); dump_top(re_p2, "Reactome_P2")

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
if (nrow(sig) == 0) {
  cat("AVISO: nenhuma via significativa (FDR<0.05). Veja results/*_top_pvalue.csv para inspecao.\n")
  cls <- data.frame(source = character(), ID = character(), Description = character(),
                    p1 = numeric(), p2 = numeric(), in_p1 = logical(), in_p2 = logical(),
                    class = factor(levels = c("High-Confidence","Hypothesis-Driven","Agnostic-Specific")))
} else {
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
}
write.table(cls, "results/pathway_convergence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
if (nrow(cls) > 0) cat("High-Confidence:", sum(cls$class == "High-Confidence"),
    "| Hypothesis-Driven:", sum(cls$class == "Hypothesis-Driven"),
    "| Agnostic-Specific:", sum(cls$class == "Agnostic-Specific"), "\n")

cat("\n=== Pronto: tabelas em results/ (figuras em plot_pathways.R) ===\n")
