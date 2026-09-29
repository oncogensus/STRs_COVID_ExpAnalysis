############################################################
# GSE188847 - DIFFERENTIAL EXPRESSION WITH edgeR
############################################################

library(edgeR)

cat("\n")
cat("============================================================\n")
cat("GSE188847 - DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("============================================================\n\n")

############################################################
# 1. FILES
############################################################

matrix_file <- "GSE188847_gene_counts.tsv"
metadata_file <- "GSE188847_gene_metadata.tsv"

############################################################
# 2. READ MATRIX
############################################################

cat("Reading matrix...\n")

counts <- read.delim(
    matrix_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

############################################################
# 3. SPLIT IDS AND COUNTS MATRIX
############################################################

gene_id <- counts$gene_id
gene_name <- counts$gene_name

count_matrix <- counts[, -(1:2)]

rownames(count_matrix) <- gene_id

count_matrix <- as.matrix(count_matrix)

storage.mode(count_matrix) <- "numeric"

############################################################
# 4. METADATA
############################################################

metadata <- read.delim(
    metadata_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample

############################################################
# Check the correspondence
############################################################

if (!all(colnames(count_matrix) == rownames(metadata))) {

    cat("\nERROR: sample order does not match!\n")

    print(colnames(count_matrix))
    print(rownames(metadata))

    stop("Fix the correspondence between matrix and metadata.")
}

############################################################
# 5. GROUPS
############################################################

metadata$group <- factor(
    metadata$group,
    levels = c(
        "CONTROL",
        "COVID",
        "ICUVENT"
    )
)

cat("\n")
cat("Samples per group:\n")
print(table(metadata$group))

############################################################
# 6. DGEList
############################################################

y <- DGEList(
    counts = count_matrix,
    group = metadata$group
)

############################################################
# 7. FILTERING
############################################################

cat("\n")
cat("Genes before filtering:", nrow(y), "\n")

keep <- filterByExpr(
    y,
    group = metadata$group
)

y <- y[keep, , keep.lib.sizes = FALSE]

cat(
    "Genes after filtering:",
    nrow(y),
    "\n"
)

############################################################
# 8. NORMALIZATION
############################################################

y <- calcNormFactors(
    y,
    method = "TMM"
)

############################################################
# 9. FUNCTION TO PERFORM DEG
############################################################

fazer_DEG <- function(
    y,
    metadata,
    group1,
    group2,
    nome
) {

    cat("\n")
    cat("============================================================\n")
    cat("COMPARISON:", group1, "vs", group2, "\n")
    cat("============================================================\n")

    ########################################################
    # Select groups
    ########################################################

    keep_samples <- metadata$group %in% c(
        group1,
        group2
    )

    y_sub <- y[, keep_samples]

    group <- droplevels(
        metadata$group[keep_samples]
    )

    ########################################################
    # Define reference
    ########################################################

    group <- factor(
        group,
        levels = c(
            group2,
            group1
        )
    )

    y_sub$samples$group <- group

    ########################################################
    # Design
    ########################################################

    design <- model.matrix(
        ~ group
    )

    ########################################################
    # Estimate dispersion
    ########################################################

    y_sub <- estimateDisp(
        y_sub,
        design
    )

    ########################################################
    # GLM
    ########################################################

    fit <- glmQLFit(
        y_sub,
        design
    )

    ########################################################
    # TEST
    ########################################################

    qlf <- glmQLFTest(
        fit,
        coef = 2
    )

    ########################################################
    # RESULTS
    ########################################################

    res <- topTags(
        qlf,
        n = Inf
    )$table

    ########################################################
    # Add gene information
    ########################################################

    res$gene_id <- rownames(res)

    res$gene_name <- gene_name[
        match(
            res$gene_id,
            gene_id
        )
    ]

    ########################################################
    # Direction
    ########################################################

    res$Significant <- ifelse(
        res$FDR < 0.05 &
        abs(res$logFC) > 1,
        "Yes",
        "No"
    )

    res$Direction <- "Not_DE"

    res$Direction[
        res$Significant == "Yes" &
        res$logFC > 1
    ] <- paste0(
        "Up_",
        group1
    )

    res$Direction[
        res$Significant == "Yes" &
        res$logFC < -1
    ] <- paste0(
        "Up_",
        group2
    )

    ########################################################
    # Reorder
    ########################################################

    res <- res[
        ,
        c(
            "gene_id",
            "gene_name",
            "logFC",
            "logCPM",
            "F",
            "PValue",
            "FDR",
            "Significant",
            "Direction"
        )
    ]

    ########################################################
    # Create directory
    ########################################################

    dir.create(
        "DEG_results",
        showWarnings = FALSE
    )

    ########################################################
    # Save COMPLETE result
    ########################################################

    write.table(
        res,
        file = paste0(
            "DEG_results/DEG_",
            nome,
            "_completo.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # NON-ADJUSTED result
    ########################################################

    res_p <- res[
        res$PValue < 0.05 &
        abs(res$logFC) > 1,
        ]

    write.table(
        res_p,
        file = paste0(
            "DEG_results/DEG_",
            nome,
            "_Pvalue.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # FDR result
    ########################################################

    res_fdr <- res[
        res$FDR < 0.05 &
        abs(res$logFC) > 1,
        ]

    write.table(
        res_fdr,
        file = paste0(
            "DEG_results/DEG_",
            nome,
            "_FDR.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # UP
    ########################################################

    up <- res[
        res$FDR < 0.05 &
        res$logFC > 1,
        ]

    write.table(
        up,
        file = paste0(
            "DEG_results/DEG_",
            nome,
            "_UP.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # DOWN
    ########################################################

    down <- res[
        res$FDR < 0.05 &
        res$logFC < -1,
        ]

    write.table(
        down,
        file = paste0(
            "DEG_results/DEG_",
            nome,
            "_DOWN.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # SUMMARY
    ########################################################

    cat("\nGenes tested:", nrow(res), "\n")

    cat(
        "DEGs with P < 0.05 and |logFC| > 1:",
        nrow(res_p),
        "\n"
    )

    cat(
        "DEGs with FDR < 0.05 and |logFC| > 1:",
        nrow(res_fdr),
        "\n"
    )

    cat(
        "UP:",
        nrow(up),
        "\n"
    )

    cat(
        "DOWN:",
        nrow(down),
        "\n"
    )

    cat("\nTop genes:\n")

    print(
        head(
            res_fdr,
            10
        )
    )

    return(res)
}

############################################################
# 10. COMPARISONS
############################################################

res_COVID_CONTROL <- fazer_DEG(
    y,
    metadata,
    "COVID",
    "CONTROL",
    "COVID_vs_CONTROL"
)

res_ICUVENT_CONTROL <- fazer_DEG(
    y,
    metadata,
    "ICUVENT",
    "CONTROL",
    "ICUVENT_vs_CONTROL"
)

res_ICUVENT_COVID <- fazer_DEG(
    y,
    metadata,
    "ICUVENT",
    "COVID",
    "ICUVENT_vs_COVID"
)

############################################################
# FINAL
############################################################

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\nFiles saved in:\n")
cat("DEG_results/\n\n")

system(
    "ls -lh DEG_results/"
)

cat("\n")