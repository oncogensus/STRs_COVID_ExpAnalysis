#!/usr/bin/env Rscript

############################################################
# EXPRESSÃO DIFERENCIAL - GSE157103
#
# Comparações:
#   1. COVID:    ICU vs NonICU
#   2. NONCOVID: ICU vs NonICU
#
# Entrada:
#   Expected Counts (RSEM)
#
# Método:
#   limma + voom
#
# Critérios de DEG:
#   FDR < 0.05
#   |log2FC| > 1
############################################################

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

############################################################
# CONFIGURAÇÕES
############################################################

dir.create("DEG_orientador", showWarnings = FALSE)

############################################################
# FUNÇÃO PRINCIPAL
############################################################

fazer_DEG <- function(
    matriz_file,
    metadata_file,
    nome
) {

    cat("\n====================================================\n")
    cat("ANÁLISE:", nome, "\n")
    cat("====================================================\n")

    ########################################################
    # DIRETÓRIO
    ########################################################

    outdir <- file.path("DEG_orientador", nome)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    ########################################################
    # LER MATRIZ
    ########################################################

    cat("\nLendo matriz...\n")

    con <- gzfile(matriz_file, "rt")
    expr <- read.delim(
        con,
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
    close(con)

    ########################################################
    # LER METADATA
    ########################################################

    cat("Lendo metadata...\n")

    meta <- read.delim(
        metadata_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    ########################################################
    # CONFERIR AMOSTRAS
    ########################################################

    if (!all(colnames(expr) == meta$GSM)) {

        cat("\nATENÇÃO: ordem das amostras diferente!\n")

        meta <- meta[
            match(colnames(expr), meta$GSM),
            ,
            drop = FALSE
        ]

    }

    if (!all(colnames(expr) == meta$GSM)) {
        stop("As amostras da matriz e metadata não correspondem.")
    }

    cat("\nAmostras:", ncol(expr), "\n")
    cat("Genes:", nrow(expr), "\n")

    ########################################################
    # GRUPOS
    ########################################################

    meta$ICU_status <- factor(
        meta$ICU_status,
        levels = c("NonICU", "ICU")
    )

    cat("\nDistribuição dos grupos:\n")
    print(table(meta$ICU_status))

    ########################################################
    # SALVAR CÓPIA DA MATRIZ E METADATA
    ########################################################

    write.table(
        expr,
        file = gzfile(
            file.path(outdir, "matriz_expressao.tsv.gz")
        ),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    write.table(
        meta,
        file = file.path(outdir, "metadata.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # DGEList
    ########################################################

    cat("\nCriando objeto DGEList...\n")

    y <- DGEList(
        counts = expr
    )

    ########################################################
    # FILTRAGEM
    ########################################################

    cat("Genes antes da filtragem:", nrow(y), "\n")

    keep <- filterByExpr(
        y,
        group = meta$ICU_status
    )

    y <- y[keep, , keep.lib.sizes = FALSE]

    cat("Genes após filtragem:", nrow(y), "\n")

    ########################################################
    # NORMALIZAÇÃO
    ########################################################

    y <- calcNormFactors(y)

    ########################################################
    # DESIGN
    ########################################################

    design <- model.matrix(
        ~ ICU_status,
        data = meta
    )

    colnames(design) <- make.names(colnames(design))

    cat("\nDesign matrix:\n")
    print(design)

    ########################################################
    # VOOM
    ########################################################

    cat("\nExecutando voom...\n")

    v <- voom(
        y,
        design,
        plot = FALSE
    )

    ########################################################
    # LIMMA
    ########################################################

    fit <- lmFit(
        v,
        design
    )

    fit <- eBayes(fit)

    ########################################################
    # RESULTADOS
    ########################################################

    # O coeficiente é ICU_statusICU
    res <- topTable(
        fit,
        coef = "ICU_statusICU",
        number = Inf,
        sort.by = "P"
    )

    ########################################################
    # ADICIONAR NOME DO GENE
    ########################################################

    res$gene <- rownames(res)

    res <- res[
        ,
        c(
            "gene",
            "logFC",
            "AveExpr",
            "t",
            "P.Value",
            "adj.P.Val",
            "B"
        )
    ]

    ########################################################
    # CLASSIFICAÇÃO
    ########################################################

    res$significance <- "Not significant"

    res$significance[
        res$adj.P.Val < 0.05 &
        res$logFC > 1
    ] <- "Up"

    res$significance[
        res$adj.P.Val < 0.05 &
        res$logFC < -1
    ] <- "Down"

    ########################################################
    # SALVAR RESULTADO COMPLETO
    ########################################################

    write.table(
        res,
        file = file.path(
            outdir,
            paste0(
                "DEG_",
                nome,
                "_completo.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # DEGS SIGNIFICATIVOS
    ########################################################

    degs <- res[
        res$adj.P.Val < 0.05 &
        abs(res$logFC) > 1,
        ,
        drop = FALSE
    ]

    write.table(
        degs,
        file = file.path(
            outdir,
            paste0(
                "DEGs_",
                nome,
                "_FDR0.05_log2FC1.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    ########################################################
    # RESUMO
    ########################################################

    n_up <- sum(
        res$adj.P.Val < 0.05 &
        res$logFC > 1
    )

    n_down <- sum(
        res$adj.P.Val < 0.05 &
        res$logFC < -1
    )

    cat("\n----------------------------------------------------\n")
    cat("RESULTADO\n")
    cat("----------------------------------------------------\n")

    cat("Genes testados:", nrow(res), "\n")
    cat("DEGs totais:", nrow(degs), "\n")
    cat("Up em ICU:", n_up, "\n")
    cat("Down em ICU:", n_down, "\n")

    ########################################################
    # PCA
    ########################################################

    pdf(
        file.path(outdir, paste0("PCA_", nome, ".pdf")),
        width = 8,
        height = 6
    )

    plotMDS(
        y,
        labels = meta$ICU_status,
        col = as.numeric(meta$ICU_status),
        main = paste("MDS -", nome)
    )

    dev.off()

    ########################################################
    # MA PLOT
    ########################################################

    pdf(
        file.path(outdir, paste0("MA_", nome, ".pdf")),
        width = 8,
        height = 6
    )

    plotMD(
        fit,
        column = 1,
        status = res$significance,
        main = paste("MA plot -", nome)
    )

    abline(
        h = c(-1, 1),
        lty = 2
    )

    dev.off()

    ########################################################
    # VOLCANO
    ########################################################

    pdf(
        file.path(outdir, paste0("Volcano_", nome, ".pdf")),
        width = 8,
        height = 6
    )

    plot(
        res$logFC,
        -log10(res$adj.P.Val),
        pch = 16,
        cex = 0.5,
        xlab = "log2 Fold Change",
        ylab = "-log10(FDR)",
        main = paste("Volcano plot -", nome)
    )

    abline(
        v = c(-1, 1),
        lty = 2
    )

    abline(
        h = -log10(0.05),
        lty = 2
    )

    dev.off()

    ########################################################
    # TOP 20
    ########################################################

    top20 <- head(degs, 20)

    write.table(
        top20,
        file = file.path(
            outdir,
            paste0(
                "TOP20_DEGs_",
                nome,
                ".tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\nArquivos salvos em:", outdir, "\n")

    invisible(res)
}

############################################################
# COVID
############################################################

res_COVID <- fazer_DEG(
    matriz_file =
        "COVID_ICU_NonICU_matrix.tsv.gz",

    metadata_file =
        "COVID_ICU_NonICU_metadata.tsv",

    nome =
        "COVID_ICU_vs_NonICU"
)

############################################################
# NONCOVID
############################################################

res_NONCOVID <- fazer_DEG(
    matriz_file =
        "NONCOVID_ICU_NonICU_matrix.tsv.gz",

    metadata_file =
        "NONCOVID_ICU_NonICU_metadata.tsv",

    nome =
        "NONCOVID_ICU_vs_NonICU"
)

############################################################
# FINAL
############################################################

cat("\n\n====================================================\n")
cat("ANÁLISES FINALIZADAS\n")
cat("====================================================\n")

cat("\nCOVID:\n")
cat("  50 ICU vs 50 NonICU\n")

cat("\nNONCOVID:\n")
cat("  16 ICU vs 10 NonICU\n")

cat("\nResultados em:\n")
cat("  DEG_orientador/\n")
