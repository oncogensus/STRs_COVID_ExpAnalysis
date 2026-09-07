# ============================================================================
# 2_skat_o_per_pathway.R -- SKAT-O por via KEGG/Reactome
# ----------------------------------------------------------------------------
# Para cada via (KEGG e Reactome), reune as STRs dos genes da via presentes na
# matrix de dosagem e roda SKAT-O (burden + SKAT, rho otimo). Relata os 3
# p-values por via, rho e correcao BH dentro de cada base de vias.
#
# Fontes de vias:
#   --gene-set-source clusterprofiler (default) | msigdbr
#     * clusterprofiler: clusterProfiler::download_KEGG("hsa") / download_Reactome("hsa")
#     * msigdbr:  msigdbr C2 CP:KEGG / CP:REACTOME (fallback robusto, sem internet)
#   O mapa via x gene e cacheado em --cache-dir (reutilizado em rodadas).
#
# Mapeamento de IDs: converte o ID da via (Entrez) para o espaco de gene_id do
# str_meta (autodetectado: ensembl/entrez/symbol) via org.Hs.eg.db.
#
# Uso:
#   Rscript 2_skat_o_per_pathway.R [--strategy gwas_burden|rna_burden|full]
#                                  [--background <arquivo>] [--out-dir <dir>]
#                                  [--gene-set-source clusterprofiler|msigdbr]
#                                  [--cache-dir <dir>] [--min-strs <n>]
# ============================================================================
source("00_common.R")

cmd_args <- commandArgs(trailingOnly = TRUE)
strategy       <- get_opt(cmd_args, "--strategy", "gwas_burden")
background_opt <- get_opt(cmd_args, "--background")
out_dir_opt    <- get_opt(cmd_args, "--out-dir")
gs_source      <- get_opt(cmd_args, "--gene-set-source", "clusterprofiler")
cache_dir_opt  <- get_opt(cmd_args, "--cache-dir")
min_strs       <- as.integer(get_opt(cmd_args, "--min-strs", "2"))
kernel         <- get_opt(cmd_args, "--kernel", "linear.weighted")

if (strategy == "gwas_burden") {
  out_dir <- file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_skat_o_pathways")
} else {
  out_dir <- file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_skat_o_pathways")
}
if (!is.null(out_dir_opt)) out_dir <- out_dir_opt
cache_dir <- ifelse(!is.null(cache_dir_opt), cache_dir_opt,
                    file.path(out_dir, "pathway_cache"))

cat("=== 2_skat_o_per_pathway [strategy:", strategy,
    "] [gene-set-source:", gs_source, "] ===\n")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Deteccao do tipo de gene_id (ensembl/entrez/symbol)
# ---------------------------------------------------------------------------
guess_id_type <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (!length(ids)) return("unknown")
  ens <- sum(grepl("^ENSG", ids)) / length(ids)
  if (ens > 0.5) return("ensembl")
  if (all(grepl("^[0-9]+$", ids))) return("entrez")
  "symbol"
}

# ---------------------------------------------------------------------------
# 2. Carrega mapa via x gene (com cache em TSV)
# ---------------------------------------------------------------------------
ensembl_db_available <- function() {
  requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)
}

entrez_to_geneid <- function(entrez) {
  entrez <- as.character(unique(entrez))
  entrez <- entrez[grepl("^[0-9]+$", entrez)]
  if (!length(entrez)) return(data.frame(entrez = character(0), gene_id = character(0)))
  if (!ensembl_db_available())
    stop("org.Hs.eg.db/AnnotationDbi ausentes. Instale via BiocManager::install(c('org.Hs.eg.db','AnnotationDbi'))")
  suppressMessages({
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez,
                                 columns = "ENSEMBL", keytype = "ENTREZID")
  })
  map <- map[!is.na(map$ENSEMBL), , drop = FALSE]
  data.frame(entrez = as.character(map$ENTREZID), gene_id = as.character(map$ENSEMBL))
}

load_pathway_map_msigdbr <- function(subcat) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("msigdbr ausente. Instale: install.packages('msigdbr')")
  gs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                         subcategory = subcat)
  data.frame(pathway = paste0("MSIGDB_", gs$gs_name),
             pathway_name = gs$gs_name,
             gene_id = as.character(gs$ensembl_gene),
             stringsAsFactors = FALSE)
}

load_pathway_map_kegg <- function() {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("clusterProfiler ausente. Instale: BiocManager::install('clusterProfiler')")
  if (!ensembl_db_available())
    stop("org.Hs.eg.db/AnnotationDbi ausentes para mapear IDs da via.")
  cat("Baixando vias KEGG (hsa) via clusterProfiler...\n")
  kegg <- clusterProfiler::download_KEGG("hsa")
  # Estrutura esperada: list(KEGGPATHID2EXTID, KEGGPATHID2NAME)
  extid <- kegg$KEGGPATHID2EXTID
  if (is.null(extid)) extid <- kegg$KEGGPATH2EXTID
  nm    <- kegg$KEGGPATHID2NAME
  if (is.null(nm)) nm <- kegg$KEGGPATH2NAME
  if (is.null(extid) || is.null(nm))
    stop("Formato inesperado de download_KEGG. Use 'msigdbr' como gene-set-source.")
  colnames(extid) <- c("pathway", "extid")
  colnames(nm)    <- c("pathway", "name")
  ent <- entrez_to_geneid(sub("^[^:]+:", "", as.character(extid$extid)))
  extid$gene_id <- ent$gene_id[match(as.character(extid$extid), ent$entrez)]
  extid <- extid[!is.na(extid$gene_id), , drop = FALSE]
  extid$pathway_name <- nm$name[match(extid$pathway, nm$pathway)]
  data.frame(pathway = paste0("KEGG_", extid$pathway),
             pathway_name = extid$pathway_name,
             gene_id = as.character(extid$gene_id),
             stringsAsFactors = FALSE)
}

load_pathway_map_reactome <- function() {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("clusterProfiler ausente. Instale: BiocManager::install('clusterProfiler')")
  if (!ensembl_db_available())
    stop("org.Hs.eg.db/AnnotationDbi ausentes para mapear IDs da via.")
  cat("Baixando vias Reactome (hsa) via clusterProfiler...\n")
  rt <- clusterProfiler::download_Reactome("hsa")
  extid <- rt$REACTOMEPATHID2EXTID
  if (is.null(extid)) extid <- rt$REACTOMEPATH2EXTID
  nm <- rt$REACTOMEPATHID2NAME
  if (is.null(nm)) nm <- rt$REACTOMEPATH2NAME
  if (is.null(extid) || is.null(nm))
    stop("Formato inesperado de download_Reactome. Use 'msigdbr' como gene-set-source.")
  colnames(extid) <- c("pathway", "extid")
  colnames(nm)    <- c("pathway", "name")
  ent <- entrez_to_geneid(sub("^[^:]+:", "", as.character(extid$extid)))
  extid$gene_id <- ent$gene_id[match(as.character(extid$extid), ent$entrez)]
  extid <- extid[!is.na(extid$gene_id), , drop = FALSE]
  extid$pathway_name <- nm$name[match(extid$pathway, nm$pathway)]
  data.frame(pathway = paste0("REACTOME_", extid$pathway),
             pathway_name = extid$pathway_name,
             gene_id = as.character(extid$gene_id),
             stringsAsFactors = FALSE)
}

load_pathway_map <- function(source, tag, cache_file) {
  if (file.exists(cache_file)) {
    cat("Cache de vias encontrado:", cache_file, "\n")
    return(fread(cache_file, header = TRUE, sep = "\t",
                 data.table = FALSE)[, c("pathway", "pathway_name", "gene_id")])
  }
  m <- switch(source,
    clusterprofiler = if (tag == "kegg") load_pathway_map_kegg()
                     else load_pathway_map_reactome(),
    msigdbr = if (tag == "kegg") load_pathway_map_msigdbr("CP:KEGG")
              else load_pathway_map_msigdbr("CP:REACTOME"),
    stop("gene-set-source invalido: ", source))
  fwrite(as.data.frame(m), cache_file, sep = "\t")
  cat("Cached:", cache_file, "linhas:", nrow(m), "\n")
  m
}

# ---------------------------------------------------------------------------
# 3. Carrega dados e roda SKAT.O por via
# ---------------------------------------------------------------------------
inp <- load_inputs(strategy, background_opt)
obj <- build_null_model(inp$covar, inp$M_dosage)

str_genes <- unique(inp$str_meta$gene_id[!is.na(inp$str_meta$gene_id)])
id_type <- guess_id_type(str_genes)
cat("Tipo de gene_id detectado:", id_type, "| genes no str_meta:", length(str_genes), "\n")

gene_to_strs <- function(genes) {
  genes <- as.character(genes)
  unlist(lapply(genes, function(g) inp$str_meta$STRs_ID[inp$str_meta$gene_id == g]))
}

fit_pathway_table <- function(map) {
  mapa <- split(map$gene_id, map$pathway)
  path_names <- map$pathway_name[match(names(mapa), map$pathway)]
  names(path_names) <- names(mapa)
  cat("Testando", length(mapa), "vias...\n")
  res <- lapply(names(mapa), function(pw) {
    genes <- unique(mapa[[pw]])
    n_genes_all <- length(genes)
    strs <- gene_to_strs(genes)
    strs <- unique(strs[strs %in% colnames(inp$M_dosage)])
    r <- run_skat_o(strs, inp$M_dosage, obj, min_var = min_strs, kernel = kernel)
    data.frame(pathway = pw, pathway_name = path_names[[pw]],
               n_genes = n_genes_all, n_strs = length(strs), r)
  })
  out <- bind_rows(res)
  out$q_value_skat_o   <- p.adjust(out$p_value_skat_o, method = "BH")
  out$q_value_skat     <- p.adjust(out$p_value_skat, method = "BH")
  out$q_value_burden   <- p.adjust(out$p_value_burden, method = "BH")
  out <- out[order(out$p_value_skat_o), ]
  if (anyNA(out$p_value_skat_o)) out <- out[order(!is.na(out$p_value_skat_o),
                                                  out$p_value_skat_o), ]
  out
}

cache_kegg <- file.path(cache_dir, "kegg_pathway_map.tsv")
cache_rt   <- file.path(cache_dir, "reactome_pathway_map.tsv")

kegg_map <- load_pathway_map(gs_source, "kegg", cache_kegg)
rt_map   <- load_pathway_map(gs_source, "reactome", cache_rt)

kegg_res <- fit_pathway_table(kegg_map)
write_hits(kegg_res, "p_value_skat_o", "q_value_skat_o",
           "skat_o_pathways_kegg", out_dir)

rt_res <- fit_pathway_table(rt_map)
write_hits(rt_res, "p_value_skat_o", "q_value_skat_o",
           "skat_o_pathways_reactome", out_dir)

all_res <- rbind(cbind(base = "KEGG", kegg_res),
                 cbind(base = "REACTOME", rt_res))
fwrite(all_res, file.path(out_dir, "skat_o_pathways_all.tsv"), sep = "\t")
fwrite(kegg_map, file.path(out_dir, "pathway_gene_kegg.tsv"), sep = "\t")
fwrite(rt_map, file.path(out_dir, "pathway_gene_reactome.tsv"), sep = "\t")

cat("\n=== SKAT-O por via KEGG (top 15) ===\n")
print(head(kegg_res, 15))
cat("\n=== SKAT-O por via Reactome (top 15) ===\n")
print(head(rt_res, 15))
unc <- all_res[!is.na(all_res$q_value_skat_o) & all_res$q_value_skat_o < 0.05, ]
cat("\n=== SKAT-O vias corrigidas (q<0.05):", nrow(unc), "===\n")
print(unc)
cat("\n=== FIM 2_skat_o_per_pathway.R ===\n")