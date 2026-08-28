#!/usr/bin/env Rscript
## gene_crossvalidation.R  (ETAPA DE GENES — cross-validation agnostic x GWAS-filtered)
## Cruza os genes com STR-outlier de duas estrategias:
##   P1 (GWAS-filtrado): outliers DBSCAN restritos a genes COVID GWAS-suggestivos (p<1e-5)
##   P2 (agnostico)    : outliers DBSCAN genome-wide
## Gera results/gene_convergence.tsv/.csv (classes) + fig_gene_upset.png (UpSet).
suppressMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(org.Hs.eg.db)
  library(UpSetR)
})

ROOT <- normalizePath("../../..")
p1_file <- "../dbscan_subset/results/covid_suggestive_genes_with_outlier_STRs.tsv"
p2_file <- "../litcovid_validation/data/outlier_genes.txt"
catalog <- file.path(ROOT, "samples/STRs_analysis_dataset.tsv")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

ORIGINAL_8 <- c("ROBO2","ANK3","CDH12","NKAIN2","SEMA6D","KCNH1","KCNQ5","ST6GALNAC3")

## ======================================================================
## 1. LISTAS DE GENES (SIMBOLOS)
## ======================================================================
p1 <- read_delim(p1_file, delim = "\t")
p1_sym <- unique(toupper(str_trim(p1$gene)))
cat("[P1] genes sugestivos com STR-outlier (GWAS-filtrado):", length(p1_sym), "\n")

p2 <- read_delim(p2_file, delim = "\t")
p2_sym <- unique(toupper(str_trim(p2$gene_name)))
cat("[P2] genes genome-wide com STR-outlier (agnostico):", length(p2_sym), "\n")

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
cat("Genes cross-validados (P1 & P2):", length(shared), "\n")

gene_conv <- data.frame(
  gene  = all_sym,
  entrez = as.character(entrez_map[match(all_sym, names(entrez_map))]),
  in_p1 = all_sym %in% p1_sym,
  in_p2 = all_sym %in% p2_sym,
  stringsAsFactors = FALSE) %>% mutate(
  strategy = case_when(in_p1 & in_p2 ~ "Cross-validated",
                       in_p1            ~ "GWAS-filtered-only",
                       in_p2            ~ "Agnostic-only",
                       TRUE             ~ "None"),
  original_8 = gene %in% ORIGINAL_8)

write.table(gene_conv, "results/gene_convergence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(gene_conv, "results/gene_convergence.csv", row.names = FALSE)

cat("=== Gene-level cross-validation: classes ===\n")
cat("Cross-validated   :", sum(gene_conv$strategy == "Cross-validated"), "\n")
cat("GWAS-filtered-only:", sum(gene_conv$strategy == "GWAS-filtered-only"), "\n")
cat("Agnostic-only     :", sum(gene_conv$strategy == "Agnostic-only"), "\n")
cat("Dentre os 8 genes STR originais presentes:", sum(gene_conv$original_8), "\n")

## ======================================================================
## 4. FIGURA UpSet (sobreposicao P1 x P2)
## ======================================================================
lst <- list(`GWAS-filtered` = p1_sym, `Agnostic` = p2_sym)
png("results/fig_gene_upset.png", width = 8, height = 6, units = "in", res = 300, bg = "white")
upset(fromList(lst), sets = c("GWAS-filtered", "Agnostic"),
      order.by = "freq", number.angles = 30,
      text.scale = c(1.4, 1.2, 1.2, 1, 1.4, 1.2),
      mainbar.y.label = "Genes por intersecao", sets.x.label = "Total por estrategia")
dev.off()
cat("Fig gene UpSet ok: results/fig_gene_upset.png\n")

cat("\n=== Pronto: results/gene_convergence.tsv + fig_gene_upset.png ===\n")
