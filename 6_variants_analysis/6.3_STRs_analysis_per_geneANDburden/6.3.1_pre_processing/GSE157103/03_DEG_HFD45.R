# GSE157103 - Differential expression analysis using HFD-45
#
# Analysis:
#   HFD-45 = 0 vs HFD-45 > 0
#   Model 1: HFD group
#   Model 2: HFD group adjusted for ICU status
#
# DEG criteria:
#   FDR < 0.05
#   |log2FC| > 1
#

library(edgeR)

# 1. FILES

matrix_file <- "COVID_HFD45_matrix.tsv.gz"
metadata_file <- "COVID_HFD45_metadata.tsv"

# 2. READ COUNT MATRIX

counts <- read.delim(
  gzfile(matrix_file),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

counts <- as.matrix(counts)
mode(counts) <- "numeric"

# 3. READ METADATA

meta <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 4. MATCH SAMPLE ORDER

common_samples <- intersect(colnames(counts), meta$GSM)

counts <- counts[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$GSM), , drop = FALSE]

rownames(meta) <- meta$GSM

stopifnot(identical(colnames(counts), rownames(meta)))


# 5. CLINICAL VARIABLES

meta$HFD_group <- factor(
  meta$HFD_group,
  levels = c("HFD_positive", "HFD0")
)

meta$ICU_status <- factor(
  meta$ICU_status,
  levels = c("no", "yes")
)

# 6. CREATE DGEList

y <- DGEList(
  counts = counts
)

# 7. FILTER LOWLY EXPRESSED GENES

keep <- filterByExpr(
  y,
  group = meta$HFD_group
)

y <- y[keep, , keep.lib.sizes = FALSE]

# 8. TMM NORMALIZATION

y <- calcNormFactors(
  y,
  method = "TMM"
)

# 9. MODEL 1 - HFD-45

design_hfd <- model.matrix(
  ~ HFD_group,
  data = meta
)

y <- estimateDisp(
  y,
  design_hfd
)

fit_hfd <- glmQLFit(
  y,
  design_hfd,
  robust = TRUE
)

qlf_hfd <- glmQLFTest(
  fit_hfd,
  coef = "HFD_groupHFD0"
)

res_hfd <- topTags(
  qlf_hfd,
  n = Inf
)$table

res_hfd$Gene <- rownames(res_hfd)

res_hfd <- res_hfd[, c(
  "Gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)]

# 10. CLASSIFY DEGs - HFD-45

res_hfd$Significant <- ifelse(
  res_hfd$FDR < 0.05 &
    abs(res_hfd$logFC) > 1,
  "Yes",
  "No"
)

res_hfd$Direction <- ifelse(
  res_hfd$Significant == "No",
  "Not_DE",
  ifelse(
    res_hfd$logFC > 0,
    "Up_HFD0",
    "Down_HFD0"
  )
)

# 11. SAVE HFD-45 RESULTS

write.table(
  res_hfd,
  file = "DEG_HFD45_bruto.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_hfd[res_hfd$Significant == "Yes", ],
  file = "DEG_HFD45_bruto_significativos.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 12. MODEL 2 - ADJUSTED FOR ICU

design_adj <- model.matrix(
  ~ ICU_status + HFD_group,
  data = meta
)

y <- estimateDisp(
  y,
  design_adj
)

fit_adj <- glmQLFit(
  y,
  design_adj,
  robust = TRUE
)

qlf_adj <- glmQLFTest(
  fit_adj,
  coef = "HFD_groupHFD0"
)

res_adj <- topTags(
  qlf_adj,
  n = Inf
)$table

res_adj$Gene <- rownames(res_adj)

res_adj <- res_adj[, c(
  "Gene",
  "logFC",
  "logCPM",
  "F",
  "PValue",
  "FDR"
)]

# 13. CLASSIFY DEGs - ADJUSTED MODEL

res_adj$Significant <- ifelse(
  res_adj$FDR < 0.05 &
    abs(res_adj$logFC) > 1,
  "Yes",
  "No"
)

res_adj$Direction <- ifelse(
  res_adj$Significant == "No",
  "Not_DE",
  ifelse(
    res_adj$logFC > 0,
    "Up_HFD0",
    "Down_HFD0"
  )
)

# 14. SAVE ADJUSTED RESULTS

write.table(
  res_adj,
  file = "DEG_HFD45_ajustado_ICU.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_adj[res_adj$Significant == "Yes", ],
  file = "DEG_HFD45_ajustado_ICU_significativos.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 15. SAVE EDGE R OBJECTS

save(
  y,
  meta,
  design_hfd,
  fit_hfd,
  res_hfd,
  design_adj,
  fit_adj,
  res_adj,
  file = "DEG_HFD45_edgeR_objects.RData"
)

# 16. SUMMARY

cat("\n")
cat("============================================\n")
cat("GSE157103 - HFD-45 differential expression\n")
cat("============================================\n\n")

cat("Samples:", ncol(y), "\n")
cat("Genes after filtering:", nrow(y), "\n\n")

cat("HFD groups:\n")
print(table(meta$HFD_group))

cat("\nICU status:\n")
print(table(meta$ICU_status))

cat("\nHFD-45 model:\n")
cat(
  "DEGs (FDR < 0.05 and |logFC| > 1):",
  sum(res_hfd$Significant == "Yes"),
  "\n"
)

cat("\nAdjusted model:\n")
cat(
  "DEGs (FDR < 0.05 and |logFC| > 1):",
  sum(res_adj$Significant == "Yes"),
  "\n"
)

cat("\nDirection - HFD model:\n")
print(table(res_hfd$Direction))

cat("\nDirection - adjusted model:\n")
print(table(res_adj$Direction))

cat("\nAnalysis completed.\n")
