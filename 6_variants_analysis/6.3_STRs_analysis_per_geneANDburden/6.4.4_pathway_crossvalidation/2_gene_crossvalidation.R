#!/usr/bin/env Rscript
## gene_crossvalidation.R  (ETAPA DE GENES — cross-validation de duas estrategias)
## Cruza os genes com STR-outlier de duas estrategias (ex.: GWAS-filtrado x RNA).
##   P1 (GWAS-filtrado): outliers DBSCAN restritos a genes COVID GWAS-suggestivos (p<1e-5)
##   P2 (RNA)          : outliers DBSCAN em genes DEGs de RNA-seq
## Gera results/<out>_gene_convergence.tsv/.csv (classes). Figuras ficam em 3_plot_pathways.R.
##
## Uso:
##   Rscript 2_gene_crossvalidation.R [--p1-file P] [--p2-file P]
##                                   [--p1-label L] [--p2-label L] [--out NAME]
##                                   [--p1-sig-only]
## Exemplo (GWAS x RNA):
##   Rscript 2_gene_crossvalidation.R \
##     --p2-file ../6.4.2_dbscan_subset_RNA/6.4.2.2_RNA_matrix/results/rna_outlier_genes.tsv \
##     --p1-label GWAS --p2-label RNA --out gwas_rna
suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(org.Hs.eg.db)
})

ROOT <- normalizePath("../../..")
p1_file <- "../6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file <- "../data/outlier_genes.txt"

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
args <- commandArgs(trailingOnly = TRUE)
p1_file <- get_opt(args, "--p1-file", p1_file)
p2_file <- get_opt(args, "--p2-file", p2_file)
p1_lab  <- get_opt(args, "--p1-label", "P1")
p2_lab  <- get_opt(args, "--p2-label", "P2")
out_tag <- get_opt(args, "--out", "gene_convergence")
p1_sig_only <- "--p1-sig-only" %in% args

cat("[args] p1_file =", p1_file, "\n")
cat("[args] p2_file =", p2_file, "\n")
cat("[args] labels  =", p1_lab, "vs", p2_lab, "\n")
cat("[args] out tag =", out_tag, "\n")
cat("[args] p1_sig_only =", p1_sig_only, "\n")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

ORIGINAL_8 <- c("ROBO2","ANK3","CDH12","NKAIN2","SEMA6D","KCNH1","KCNQ5","ST6GALNAC3")

## ======================================================================
## 1. LISTAS DE GENES (SIMBOLOS)
## ======================================================================
read_gene_list <- function(path, sig_only = FALSE) {
  d <- read_delim(path, delim = "\t")
  if (sig_only) {
    sc <- grep("^gwas_significance$", colnames(d), value = TRUE)
    if (length(sc) == 0) stop("coluna 'gwas_significance' nao encontrada em: ", path)
    d <- d[d[[sc]] == "significant", , drop = FALSE]
    if (nrow(d) == 0) stop("nenhum gene com gwas_significance='significant' em: ", path)
  }
  gc <- grep("^gene(_name)?$", colnames(d), value = TRUE)[1]
  if (is.na(gc)) stop("coluna 'gene'/'gene_name' nao encontrada em: ", path)
  unique(toupper(str_trim(d[[gc]])))
}

p1_sym <- read_gene_list(p1_file, sig_only = p1_sig_only)
if (p1_sig_only) {
  cat(sprintf("[%s] genes significativos (p<5e-8) com STR-outlier: %d\n",
              p1_lab, length(p1_sym)))
} else {
  cat(sprintf("[%s] genes com STR-outlier: %d\n", p1_lab, length(p1_sym)))
}

p2_sym <- read_gene_list(p2_file)
cat(sprintf("[%s] genes com STR-outlier: %d\n", p2_lab, length(p2_sym)))

## ======================================================================
## 2. MAPEAMENTO SYMBOL -> ENTREZ
## ======================================================================
sym2entrez <- function(syms) {
  mapIds(org.Hs.eg.db, keys = syms, column = "ENTREZID",
         keytype = "SYMBOL", multiVals = "first")
}
all_sym <- union(p1_sym, p2_sym)
entrez_map <- sym2entrez(all_sym)
entrez_map <- entrez_map[!is.na(entrez_map)]

## ======================================================================
## 3. CLASSIFICACAO DE CROSS-VALIDATION (nivel gene)
## ======================================================================
shared <- intersect(p1_sym, p2_sym)
cat(sprintf("Genes cross-validados (%s & %s): %d\n", p1_lab, p2_lab, length(shared)))

p1_only_lab <- paste0(p1_lab, "-only")
p2_only_lab <- paste0(p2_lab, "-only")

gene_conv <- data.frame(
  gene  = all_sym,
  entrez = as.character(entrez_map[match(all_sym, names(entrez_map))]),
  in_p1 = all_sym %in% p1_sym,
  in_p2 = all_sym %in% p2_sym,
  stringsAsFactors = FALSE) %>% mutate(
  strategy = case_when(in_p1 & in_p2 ~ "Cross-validated",
                       in_p1            ~ p1_only_lab,
                       in_p2            ~ p2_only_lab,
                       TRUE             ~ "None"),
  original_8 = gene %in% ORIGINAL_8)

out_tsv <- sprintf("results/%s.tsv", out_tag)
out_csv <- sprintf("results/%s.csv", out_tag)
write.table(gene_conv, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(gene_conv, out_csv, row.names = FALSE)

cat("=== Gene-level cross-validation: classes ===\n")
cat("Cross-validated :", sum(gene_conv$strategy == "Cross-validated"), "\n")
cat(p1_only_lab, ":", sum(gene_conv$strategy == p1_only_lab), "\n")
cat(p2_only_lab, ":", sum(gene_conv$strategy == p2_only_lab), "\n")
cat("Dentre os 8 genes STR originais presentes:", sum(gene_conv$original_8), "\n")

cat(sprintf("\n=== Pronto: %s ===\n", out_tsv))
