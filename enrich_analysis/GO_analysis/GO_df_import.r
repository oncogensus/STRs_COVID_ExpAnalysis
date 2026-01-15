## Libs
library(stringr)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)



## Data import - STRs OUTLIERS (foreground with p_adj < 0.05)
df_strs_outliers <- read_delim(
  "../samples/gp_global/merged/outliers_annotated_complete.tsv",
  delim = "\t"
)
cat("Lines in df_strs_outliers:", nrow(df_strs_outliers), "\n")



## Data import - scRNA-Seq
df_genes <- read_csv("../samples/gp_global/merged/merged_scovid.csv")
cat("Lines in df_genes:", nrow(df_genes), "\n")
df_genes <- df_genes %>%
  select(-Dataset, -Tissue, -Accession)



## ==================================================================
## CLEANING
## ==================================================================


## To df_strs_outliers: limpar gene_id
df_strs_outliers <- df_strs_outliers %>%
  mutate(
    gene_id = str_trim(gene_id),
    gene_id = str_replace_all(gene_id, '^"|"$', ''),
    gene_id = toupper(gene_id)
  )



## To df_genes: clean gene_name
df_genes <- df_genes %>%
  mutate(
    `Gene Symbol` = str_trim(`Gene Symbol`),
    `Gene Symbol` = str_replace_all(`Gene Symbol`, '^"|"$|^--', ''),
    `Gene Symbol` = toupper(`Gene Symbol`)
  )



## Filter only non-empty values
df_strs_outliers <- df_strs_outliers %>%
  filter(!is.na(gene_id) & gene_id != "")

df_genes <- df_genes %>%
  filter(!is.na(`Gene Symbol`) & `Gene Symbol` != "")



cat("df_strs_outliers after cleaning:", nrow(df_strs_outliers), "\n")
cat("df_genes after cleaning:", nrow(df_genes), "\n")



## ==================================================================
## MAPPING: gene_id (ENSEMBL) → EntrezID
## ==================================================================


## Map gene_id (ENSEMBL) to Entrez ID for df_strs_outliers
df_strs_outliers <- df_strs_outliers %>%
  mutate(
    EntrezID = mapIds(
      org.Hs.eg.db,
      keys    = gene_id,
      column  = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  filter(!is.na(EntrezID))


cat("df_strs_outliers after Ensembl mapping:", nrow(df_strs_outliers), "\n")



## Map Gene Symbol to Entrez ID for df_genes (for background)
df_genes <- df_genes %>%
  mutate(
    EntrezID = mapIds(
      org.Hs.eg.db,
      keys    = `Gene Symbol`,
      column  = "ENTREZID",
      keytype = "ALIAS",
      multiVals = "first"
    )
  ) %>%
  filter(!is.na(EntrezID))


cat("df_genes after mapping:", nrow(df_genes), "\n")



## ==================================================================
## DEFINE BACKGROUND AND FOREGROUND FOR GO ENRICHMENT
## ==================================================================


## Background: ALL genes expressed in the lung (regardless of whether they have STRs)
## This is the real biological universe of the lung tissue
background_gene_list <- unique(df_genes$EntrezID)
cat("Number of background genes (ALL lung-expressed genes):", length(background_gene_list), "\n")



## Foreground: STRs with statistical significance (p_adj < 0.05)
## Filter only genes that are also in the background
case_gene_list <- unique(
  df_strs_outliers$EntrezID[
    df_strs_outliers$p_adj < 0.05 & !is.na(df_strs_outliers$p_adj)
  ]
)
## Ensure that the foreground is contained within the background
case_gene_list <- intersect(case_gene_list, background_gene_list)

cat("Number of foreground genes (p_adj < 0.05 in background):", length(case_gene_list), "\n")



## Check proportion
proportion <- length(case_gene_list) / length(background_gene_list)
cat("Proportion foreground/background:", round(proportion, 3), "\n")



## ==================================================================
## GO ENRICHMENT ANALYSIS
## ==================================================================


## GO enrichment with Entrez IDs ready
enrich_result <- enrichGO(
  gene         = case_gene_list,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL",
  universe     = background_gene_list,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 1.0
)



## Check if there are results before simplifying
if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
  
  ## Simplify terms for easier visualization
  simp_resultsGO <- clusterProfiler::simplify(enrich_result)
  
  # Create directory if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results", recursive = TRUE)
    cat("Directory 'results' created\n")
  }
  
  # Save RDS object
  saveRDS(simp_resultsGO, file = "results/simp_resultsGO.rds")
  
  ## Transform results into data frame and save CSV for inspection
  dat.results <- data.frame(simp_resultsGO, stringsAsFactors = FALSE)
  write.csv(dat.results, file = "results/GOenrich_p-value0.05_allSTRs.csv", row.names = FALSE)
  
  cat("\n✅ GO enrichment completed successfully!\n")
  cat("Number of enriched terms:", nrow(dat.results), "\n")
  
} else {
  cat("\n⚠️ WARNING: No significantly enriched GO terms found.\n")
  cat("This may indicate:\n")
  cat("  1. True negative result (no enrichment)\n")
  cat("  2. Small foreground set\n")
  cat("  3. Need to adjust pvalueCutoff or qvalueCutoff\n")
}
