# ============================================================================
# 3_enrichment_clusterprofiler.R -- enriquecimento de vias (KEGG + Reactome)
# ----------------------------------------------------------------------------
# Recebe uma lista de genes de interesse (Ensembl) e, opcionalmente, um
# universo (background). Testa enriquecimento em vias KEGG e Reactome.
#
# Engine 1 (default): clusterProfiler
#   - enrichKEGG(gene, universe, organism="hsa", keyType="ncbi-geneid")
#   - ReactomePA::enrichPathway(gene, organism="hsa") (usa Entrez)
#   Necessita de org.Hs.eg.db p/ bitr(ENSEMBL->ENTREZID) e internet no nodo.
# Engine 2 (fallback automático): teste hipergeometrico com
#   mapas de vias cacheados (mesmos TSVs do 2_skat_o_per_pathway.R),
#   em espaco Ensembl, com correcao BH.
#
# Também cruza os genes de interesse com a via KEGG COVID-19 (hsa05171) e
# com vias COVID-19 da Reactome, e reporta STRs cobertos (via norm_file).
#
# Uso:
#   Rscript 3_enrichment_clusterprofiler.R --gene-file <tsv|txt|rds> \
#       [--universe-file <arquivo>] [--label <nome>] [--out-dir <dir>] \
#       [--cache-dir <dir>] [--gene-set-source clusterprofiler|msigdbr] \
#       [--gene-id-col gene]
#
# Formato do --gene-file: TSV com coluna de gene (default 'gene') OU arquivo
# texto com um gene por linha (detectado se não houver header com colunas).
# ============================================================================
source("00_common.R")

cmd_args <- commandArgs(trailingOnly = TRUE)
gene_file     <- get_opt(cmd_args, "--gene-file", stop("--gene-file obrigatorio"))
universe_file <- get_opt(cmd_args, "--universe-file")
label         <- get_opt(cmd_args, "--label", "genes")
out_dir_opt   <- get_opt(cmd_args, "--out-dir")
cache_dir_opt <- get_opt(cmd_args, "--cache-dir")
gs_source     <- get_opt(cmd_args, "--gene-set-source", "clusterprofiler")
gene_id_col   <- get_opt(cmd_args, "--gene-id-col", "gene")

out_dir <- if (!is.null(out_dir_opt)) out_dir_opt else
  file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_enrichment")
cache_dir <- ifelse(!is.null(cache_dir_opt), cache_dir_opt,
                    file.path(out_dir, "pathway_cache"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 3_enrichment_clusterprofiler [", label, "] ===\n")

# ---------------------------------------------------------------------------
# 1. Le genes de interesse e universo (Ensembl).
# ---------------------------------------------------------------------------
read_gene_list <- function(path) {
  if (grepl("\\.rds$", path)) {
    all <- readRDS(path)
    if (is.data.frame(all)) all[[1]] else as.character(all)
  } else {
    head <- tryCatch({
      con <- file(path, open = "r")
      on.exit(close(con))
      readLines(con, n = 1)
    }, error = function(e) NULL)
    # heuristica simples: se a 1a linha tem tab e contem o nome da coluna,
    # trata como TSV; senao, um gene por linha.
    is_tsv <- !is.null(head) && grepl("\t", head) &&
      grepl(sprintf("(?i)^%s", gene_id_col), head)
    if (is_tsv) {
      fread(path, header = TRUE, sep = "\t",
            select = gene_id_col, data.table = FALSE)[[1]]
    } else {
      readLines(path, warn = FALSE)
    }
  }
}

genes <- unique(trimws(as.character(read_gene_list(gene_file))))
genes <- genes[nzchar(genes) & !is.na(genes)]
cat("Genes de interesse:", length(genes), "\n")

universe <- if (!is.null(universe_file))
  unique(trimws(as.character(read_gene_list(universe_file)))) else NULL
if (!is.null(universe)) universe <- universe[nzchar(universe)]
if (!is.null(universe) && length(universe) > 0) {
  cat("Universo (background):", length(universe), "\n")
} else {
  cat("Sem universo informado (enriquecimento sem background).\n")
  universe <- character(0)
}

# ---------------------------------------------------------------------------
# 2. Helpers de conversao ENSEMBL <-> ENTREZID (org.Hs.eg.db).
# ---------------------------------------------------------------------------
entrez_ready <- function() {
  requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)
}

ens2entrez <- function(genes_ens) {
  if (!entrez_ready()) return(NULL)
  suppressMessages({
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(genes_ens),
                          column = "ENTREZID", keytype = "ENSEMBL",
                          multiVals = "first")
  })
}

# ---------------------------------------------------------------------------
# 3. Engine clusterProfiler (enrichKEGG + ReactomePA::enrichPathway).
# ---------------------------------------------------------------------------
run_clusterprofiler <- function(genes_ens, universe_ens) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("clusterProfiler ausente")
  if (!entrez_ready()) stop("org.Hs.eg.db ausente")
  g2e <- ens2entrez(genes_ens)
  if (length(g2e) == 0) return(NULL)
  names(g2e) <- NULL
  g2e <- unique(as.character(g2e[!is.na(g2e)]))
  if (!length(g2e)) return(NULL)
  u2e <- if (length(universe_ens))
    unique(as.character(ens2entrez(universe_ens))) else NULL
  u2e <- u2e[!is.na(u2e)]
  if (!length(u2e)) u2e <- NULL

  out <- list()
  k <- tryCatch(clusterProfiler::enrichKEGG(gene = g2e, universe = u2e,
                                            organism = "hsa",
                                            keyType = "ncbi-geneid",
                                            pvalueCutoff = 1,
                                            qvalueCutoff = 1),
                error = function(e) NULL)
  if (!is.null(k) && nrow(k@result) > 0) {
    kr <- as.data.frame(k)
    out[["KEGG"]] <- data.frame(
      base = "KEGG",
      pathway = paste0(kr$ID, " - ", kr$Description),
      pathway_name = kr$Description,
      n_genes_in_pathway = kr$setSize,
      n_overlap = kr$Count,
      genes_overlap = vapply(kr$geneID, function(x)
        paste(entrez2symbol(x), collapse = ";"), character(1)),
      p_value = kr$pvalue, q_value = kr$p.adjust,
      stringsAsFactors = FALSE)
  }
  if (requireNamespace("ReactomePA", quietly = TRUE)) {
    r <- tryCatch(ReactomePA::enrichPathway(gene = g2e, organism = "hsa",
                                            universe = u2e, pvalueCutoff = 1,
                                            qvalueCutoff = 1, readable = FALSE),
                  error = function(e2) NULL)
    if (!is.null(r) && nrow(r@result) > 0) {
      kr <- as.data.frame(r)
      out[["REACTOME"]] <- data.frame(
        base = "REACTOME",
        pathway = paste0(kr$ID, " - ", kr$Description),
        pathway_name = kr$Description,
        n_genes_in_pathway = kr$setSize,
        n_overlap = kr$Count,
        genes_overlap = kr$geneID,
        p_value = kr$pvalue, q_value = kr$p.adjust,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  bind_rows(out) %>% arrange(.data$p_value)
}

entrez2symbol <- function(x) {
  if (!nzchar(x)) return(character(0))
  ids <- unlist(strsplit(as.character(x), "/"))
  suppressMessages({
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                          column = "SYMBOL", keytype = "ENTREZID",
                          multiVals = "first")
  })
}

# ---------------------------------------------------------------------------
# 4. Engine msigdbr (hipergeometrico, espaco Ensembl) -- fallback/cache.
#    get_map() e enrich_hyper() vêm de 00_common.R.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 5. Executar enriquecimento e gravar.
# ---------------------------------------------------------------------------
kegg_map <- get_map("kegg", cache_dir)
rt_map   <- get_map("reactome", cache_dir)

result <- NULL
if (gs_source == "clusterprofiler") {
  cat("Tentando enriquecimento via clusterProfiler (KEGG + Reactome)...\n")
  result <- tryCatch(run_clusterprofiler(genes, universe),
                     error = function(e) {
                       cat("Aviso: clusterProfiler falhou (", conditionMessage(e),
                           "). Usando fallback msigdbr.\n")
                       NULL
                     })
  if (is.null(result) || (is.data.frame(result) && nrow(result) == 0)) {
    cat("clusterProfiler sem resultados -- usando fallback msigdbr.\n")
    result <- NULL
  }
  if (!is.null(result)) {
    result$q_value <- if (all(is.na(result$q_value))) p.adjust(result$p_value, "BH")
                      else result$q_value
  }
}
if (is.null(result)) {
  cat("Usando mapas de vias (hipergeometrico) com os caches locais.\n")
  kegg_e   <- enrich_hyper(genes, universe, kegg_map, "KEGG")
  react_e  <- enrich_hyper(genes, universe, rt_map, "REACTOME")
  result   <- rbind(kegg_e, react_e)
  result$q_value <- p.adjust(result$p_value, method = "BH")
}

if (is.null(result) || nrow(result) == 0) {
  cat("Nenhuma via enriquecida. Verifique a lista de genes e o universo.\n")
  quit(save = "no", status = 0)
}

result <- result[order(result$p_value), ]
fwrite(result, file.path(out_dir, paste0("enrichment_", label, ".tsv")), sep = "\t")
sig_e <- result[!is.na(result$q_value) & result$q_value < 0.05, , drop = FALSE]
fwrite(sig_e, file.path(out_dir, paste0("enrichment_", label, "_significant.tsv")),
       sep = "\t")

cat("\n=== Enriquecimento (top 20) ===\n")
print(head(result, 20))
cat("\n=== Vias significativas (q<0.05):", nrow(sig_e), "===\n")
print(sig_e)

# Plot (dotplot)
if (nrow(sig_e) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plotdf <- head(sig_e[order(sig_e$q_value, sig_e$n_overlap), ], 20)
  p <- ggplot(plotdf, aes(x = reorder(pathway, -log10(q_value)),
                          y = -log10(q_value), size = n_overlap)) +
    geom_point(color = "steelblue") + coord_flip() +
    labs(title = paste("Enriquecimento de vias (BH) -", label),
         x = NULL, y = "-log10(q)") + theme_bw()
  ggsave(file.path(out_dir, paste0("dotplot_", label, ".png")), p,
         width = 9, height = max(4, 0.35 * nrow(plotdf)), dpi = 150)
}

# ---------------------------------------------------------------------------
# 6. Cruzamento explicito com vias COVID-19 (KEGG hsa05171 | Reactome).
# ---------------------------------------------------------------------------
cov_kegg_mask <- grepl("05171|COVID|CORONAVIRUS", paste(kegg_map$pathway,
                                                        kegg_map$pathway_name),
                       ignore.case = TRUE)
cov_kegg_genes <- unique(kegg_map$gene_id[cov_kegg_mask])
cov_overlap <- intersect(cov_kegg_genes, genes)
cov_rt_mask <- grepl("COVID|CORONAVIRUS", paste(rt_map$pathway, rt_map$pathway_name),
                     ignore.case = TRUE)
cov_rt_genes <- unique(rt_map$gene_id[cov_rt_mask])
cov_rt_overlap <- intersect(cov_rt_genes, genes)

# STRs cobertos (via norm_file)
cov_strs <- tryCatch({
  norm <- fread(norm_file, header = TRUE, sep = "\t", data.table = FALSE)
  cov_kegg_strs <- unique(norm$STRs_ID[norm$gene_id %in% cov_overlap])
  cov_rt_strs   <- unique(norm$STRs_ID[norm$gene_id %in% cov_rt_overlap])
  list(kegg_strs = cov_kegg_strs, rt_strs = cov_rt_strs)
}, error = function(e) list(kegg_strs = character(0), rt_strs = character(0)))

cov_tab <- data.frame(
  source = c("KEGG_hsa05171_COVID19", "Reactome_COVID19"),
  n_genes_pathway = c(length(cov_kegg_genes), length(cov_rt_genes)),
  n_genes_candidate = c(length(cov_overlap), length(cov_rt_overlap)),
  n_strs_candidate = c(length(cov_strs$kegg_strs), length(cov_strs$rt_strs)),
  genes = c(paste(cov_overlap, collapse = ";"),
            paste(cov_rt_overlap, collapse = ";")),
  strs = c(paste(cov_strs$kegg_strs, collapse = ";"),
           paste(cov_strs$rt_strs, collapse = ";")),
  stringsAsFactors = FALSE)
fwrite(cov_tab, file.path(out_dir, paste0("covid19_pathway_overlap_", label, ".tsv")),
       sep = "\t")
cat("\n=== Overlap com vias COVID-19 ===\n")
print(cov_tab)
cat("\n=== FIM 3_enrichment_clusterprofiler.R ===\n")