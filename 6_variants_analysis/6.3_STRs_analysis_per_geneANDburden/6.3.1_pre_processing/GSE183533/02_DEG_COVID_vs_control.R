#!/usr/bin/env Rscript

library(edgeR)
library(limma)

cat("\n")
cat("============================================================\n")
cat("GSE183533 - EXPRESSÃO DIFERENCIAL COM edgeR\n")
cat("============================================================\n\n")

# 1. ARQUIVOS

counts_file <- "GSE183533_gene_counts.tsv"
metadata_file <- "GSE183533_metadata.tsv"

output_dir <- "DEG_results"

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# 2. LER MATRIZ

cat("Lendo matriz...\n")

counts <- read.delim(
    counts_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# primeira coluna = gene_id
gene_ids <- counts$gene_id

counts_matrix <- counts[, -1]

rownames(counts_matrix) <- gene_ids

# converter para matriz numérica
counts_matrix <- as.matrix(counts_matrix)

storage.mode(counts_matrix) <- "numeric"

cat("Genes:", nrow(counts_matrix), "\n")
cat("Amostras:", ncol(counts_matrix), "\n\n")

# 3. LER METADATA

metadata <- read.delim(
    metadata_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# garantir mesma ordem
metadata <- metadata[
    match(colnames(counts_matrix), metadata$sample),
]

# verificar
if (!all(metadata$sample == colnames(counts_matrix))) {
    stop("ERRO: amostras da matriz e metadata não estão alinhadas.")
}

metadata$group <- factor(
    metadata$group,
    levels = c("CONTROL", "COVID")
)

cat("Amostras por grupo:\n\n")
print(table(metadata$group))

# 4. DGEList

y <- DGEList(
    counts = counts_matrix,
    group = metadata$group
)

cat("\nGenes antes da filtragem:", nrow(y), "\n")

# 5. FILTRAGEM

keep <- filterByExpr(
    y,
    group = metadata$group
)

y <- y[keep, , keep.lib.sizes = FALSE]

cat("Genes após filtragem:", nrow(y), "\n")

# 6. NORMALIZAÇÃO TMM

y <- calcNormFactors(y, method = "TMM")

cat("\nFatores de normalização TMM:\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

# 7. DESIGN

design <- model.matrix(
    ~ group,
    data = metadata
)

colnames(design) <- make.names(colnames(design))

cat("\nDesign:\n")
print(design)

# 8. ESTIMATIVA DISPERSÃO

y <- estimateDisp(y, design)

# 9. MODELO GLM

fit <- glmQLFit(
    y,
    design,
    robust = TRUE
)

# 10. CONTRASTE COVID VS CONTROL

contrast <- makeContrasts(
    groupCOVID,
    levels = design
)

qlf <- glmQLFTest(
    fit,
    contrast = contrast
)

# 11. RESULTADOS

results <- topTags(
    qlf,
    n = Inf
)$table

results$gene_id <- rownames(results)

# reorganizar
results <- results[
    ,
    c(
        "gene_id",
        "logFC",
        "logCPM",
        "F",
        "PValue",
        "FDR"
    )
]

# 12. CLASSIFICAÇÃO

results$Significant <- ifelse(
    results$FDR < 0.05 &
        abs(results$logFC) > 1,
    "Yes",
    "No"
)

results$Direction <- "Not_DE"

results$Direction[
    results$FDR < 0.05 &
        results$logFC > 1
] <- "Up_COVID"

results$Direction[
    results$FDR < 0.05 &
        results$logFC < -1
] <- "Down_COVID"

# 13. ORDENAR POR FDR

results <- results[
    order(results$FDR),
]

# 14. SEPARAR RESULTADOS

deg_fdr <- subset(
    results,
    FDR < 0.05 &
        abs(logFC) > 1
)

deg_up <- subset(
    results,
    FDR < 0.05 &
        logFC > 1
)

deg_down <- subset(
    results,
    FDR < 0.05 &
        logFC < -1
)

# 15. RESUMO

cat("\n")
cat("============================================================\n")
cat("COMPARAÇÃO: COVID vs CONTROL\n")
cat("============================================================\n\n")

cat("Genes testados:", nrow(results), "\n")

cat(
    "DEG por FDR < 0.05 e |logFC| > 1:",
    nrow(deg_fdr),
    "\n"
)

cat("UP:", nrow(deg_up), "\n")
cat("DOWN:", nrow(deg_down), "\n")

# 16. TOP GENES

cat("\nTop genes:\n\n")

print(
    head(
        deg_fdr,
        20
    )
)

# 17. SALVAR

write.table(
    results,
    file = file.path(
        output_dir,
        "DEG_COVID_vs_CONTROL_completo.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    deg_fdr,
    file = file.path(
        output_dir,
        "DEG_COVID_vs_CONTROL_FDR.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    deg_up,
    file = file.path(
        output_dir,
        "DEG_COVID_vs_CONTROL_UP.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    deg_down,
    file = file.path(
        output_dir,
        "DEG_COVID_vs_CONTROL_DOWN.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# 18. NORMALIZED COUNTS

norm_counts <- cpm(
    y,
    normalized.lib.sizes = TRUE
)

norm_counts <- data.frame(
    gene_id = rownames(norm_counts),
    norm_counts,
    check.names = FALSE
)

write.table(
    norm_counts,
    file = file.path(
        output_dir,
        "GSE183533_normalized_CPM.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# FINAL

cat("\n")
cat("============================================================\n")
cat("ANÁLISE CONCLUÍDA\n")
cat("============================================================\n\n")

cat("Arquivos salvos em:", output_dir, "\n\n")

cat(" - DEG_COVID_vs_CONTROL_completo.tsv\n")
cat(" - DEG_COVID_vs_CONTROL_FDR.tsv\n")
cat(" - DEG_COVID_vs_CONTROL_UP.tsv\n")
cat(" - DEG_COVID_vs_CONTROL_DOWN.tsv\n")
cat(" - GSE183533_normalized_CPM.tsv\n\n")
