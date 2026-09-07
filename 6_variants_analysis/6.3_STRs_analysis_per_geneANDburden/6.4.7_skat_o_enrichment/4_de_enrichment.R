# ============================================================================
# 4_de_enrichment.R -- enriquecimento de STRs em genes DE (por estudo)
# ----------------------------------------------------------------------------
# Para cada estudio RNA-seq (GSE dirs em --deg-dir com TSV/CSV de DEGs), monta
# uma tabela 2x2 (gene em DE x gene com STR no catalogo) e roda Fisher:
#   - Ha enriquecimento de genes com STRs entre os genes DE?
# Também reporta overlap DE x outlier_STRs (por gene) e, se houver mapas de
# vias cacheados, enriquecimento (hipergeometrico) das vias nos genes DE.
#
# Colunas detectadas no arquivo por estudo (flexivel):
#   Gene | gene_symbol | gene_name   (coluna de gene)
#   Significant | significance       (flag DE)
#   FDR | adj.P.Val | P.Value        (p/FDR ajustado)
#
# Uso:
#   Rscript 4_de_enrichment.R --deg-dir <dir_com_GSE> [--catalog <tsv>]
#                              [--out-dir <dir>] [--cache-dir <dir>]
#                              [--fdr-thresh 0.05] [--min-de 5]
# ============================================================================
source("00_common.R")

cat("=== 4_de_enrichment.R ===\n")

cmd_args <- commandArgs(trailingOnly = TRUE)
deg_dir       <- get_opt(cmd_args, "--deg-dir", stop("--deg-dir obrigatorio"))
catalog_opt   <- get_opt(cmd_args, "--catalog")
out_dir_opt   <- get_opt(cmd_args, "--out-dir")
cache_dir_opt <- get_opt(cmd_args, "--cache-dir")
fdr_thresh    <- as.numeric(get_opt(cmd_args, "--fdr-thresh", "0.05"))
min_de        <- as.integer(get_opt(cmd_args, "--min-de", "5"))

out_dir <- if (!is.null(out_dir_opt)) out_dir_opt else
  file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_de_enrichment")
cache_dir <- ifelse(!is.null(cache_dir_opt), cache_dir_opt,
                    file.path(out_dir, "pathway_cache"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

catalog_file <- if (!is.null(catalog_opt)) catalog_opt else
  file.path(REPO_ROOT, "samples/STRs_analysis_dataset.tsv")
if (!file.exists(catalog_file)) {
  cat("AVISO: catalogo nao encontrado em", catalog_file, "-- usando norm_file.\n")
  catalog_file <- norm_file
}

# ---------------------------------------------------------------------------
# 1. Genes com STRs (background) -- do catalogo/norm_file.
# ---------------------------------------------------------------------------
catalog <- tryCatch({
  catc <- fread(catalog_file, header = TRUE, sep = "\t", data.table = FALSE)
  catc
}, error = function(e) NULL)
if (is.null(catalog)) catalog <- data.frame()

str_genes_ens   <- if ("gene_id"   %in% colnames(catalog))
  unique(as.character(catalog$gene_id[!is.na(catalog$gene_id)])) else character(0)
str_genes_sym   <- if ("gene_name" %in% colnames(catalog))
  unique(as.character(catalog$gene_name[!is.na(catalog$gene_name)])) else character(0)
if (!length(str_genes_ens) && !length(str_genes_sym)) {
  if ("gene" %in% colnames(catalog))
    str_genes_sym <- unique(as.character(catalog$gene[!is.na(catalog$gene)]))
}
cat("Genes com STRs (background):", length(unique(c(str_genes_ens, str_genes_sym))), "\n")

# Mapa symbol -> ENSG dentro do catalogo (primeiro gene_id por gene_name)
sym_map <- if ("gene_name" %in% colnames(catalog) && length(str_genes_ens)) {
  setNames(str_genes_ens[match(str_genes_sym, str_genes_sym)],
           str_genes_sym)
} else NULL
sym_map <- unique(data.frame(gene_name = str_genes_sym,
                             gene_id   = str_genes_ens[match(str_genes_sym,
                                                             str_genes_sym)],
                             stringsAsFactors = FALSE))

# ---------------------------------------------------------------------------
# 2. Helpers de deteccao de colunas (mesma logica do cross_DEGs_STRs.py).
# ---------------------------------------------------------------------------
detect_col <- function(header, cands) {
  low <- tolower(trimws(header))
  hit <- cands[tolower(cands) %in% low]
  if (length(hit)) header[low == tolower(hit[1])][1] else NA
}

fisher_de_str <- function(de_genes, all_genes, str_genes) {
  all_genes <- unique(all_genes)
  de_genes  <- unique(de_genes)
  str_genes <- intersect(unique(str_genes), all_genes)
  de_genes  <- intersect(de_genes, all_genes)
  a <- length(intersect(de_genes, str_genes))          # DE & STR
  b <- length(setdiff(de_genes, str_genes))            # DE & !STR
  c <- length(setdiff(str_genes, de_genes))            # !DE & STR
  d <- length(setdiff(all_genes, union(de_genes, str_genes))) # !DE & !STR
  m <- matrix(c(a, b, c, d), nrow = 2,
              dimnames = list(str = c("sim", "nao"),
                              de = c("sim", "nao")))
  f <- tryCatch(fisher.test(m, alternative = "greater"), error = function(e) NULL)
  data.frame(n_de_total = length(de_genes),
             n_de_with_str = a,
             n_de_without_str = b,
             n_non_de_with_str = c,
             n_non_de_without_str = d,
             prop_de_with_str = ifelse(length(de_genes) > 0, a / length(de_genes), NA),
             prop_non_de_with_str = ifelse((c + d) > 0, c / (c + d), NA),
             odds_ratio = if (is.null(f)) NA_real_ else f$estimate,
             fisher_p = if (is.null(f)) NA_real_ else f$p.value)
}

# ---------------------------------------------------------------------------
# 3. Processa cada GSE dir.
# ---------------------------------------------------------------------------
deg_dirs <- list.dirs(deg_dir, recursive = FALSE)
if (!length(deg_dirs)) {
  # tenta tratar o proprio --deg-dir como diretorio com arquivos soltos
  deg_dirs <- deg_dir
}
if (!length(deg_dirs)) stop("Nenhum diretorio GSE encontrado em ", deg_dir)

res_fisher <- list()
res_overlap <- list()

for (d in deg_dirs) {
  gse <- basename(d)
  cat("\n--- Processando GSE:", gse, "---\n")
  files <- list.files(d, pattern = "(tsv|txt|csv)$", full.names = TRUE,
                      ignore.case = TRUE)
  if (!length(files)) {
    cat("Sem arquivos TSV/CSV em", d, "\n")
    next
  }
  f0 <- files[1]
  dt <- tryCatch({
    if (grepl("\\.csv$", f0)) fread(f0, header = TRUE, sep = ",")
    else fread(f0, header = TRUE, sep = "\t")
  }, error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) {
    cat("Falha ao ler", f0, "\n")
    next
  }
  hd <- colnames(dt)
  gene_col <- detect_col(hd, c("Gene", "gene_symbol", "gene_name", "ensembl_gene"))
  sig_col  <- detect_col(hd, c("Significant", "significance", "sig"))
  fdr_col  <- detect_col(hd, c("FDR", "adj.P.Val", "P.Value", "pvalue"))
  if (is.na(gene_col)) {
    cat("Sem coluna de gene (Gene/gene_symbol/gene_name) em", f0, "\n")
    next
  }
  dt$gene_ <- as.character(dt[[gene_col]])
  dt <- dt[!is.na(dt$gene_) & nzchar(dt$gene_), ]

  if (!is.na(sig_col)) {
    de <- dt[dt[[sig_col]] %in% c(1, TRUE, "TRUE", "yes", "Yes", "sig", "Sig", "Significant"), ]
    de_genes <- unique(dt$gene_[dt[[sig_col]] %in%
                       c(1, TRUE, "TRUE", "yes", "Yes", "sig", "Sig", "Significant")])
  } else if (!is.na(fdr_col)) {
    fv <- suppressWarnings(as.numeric(dt[[fdr_col]]))
    de_genes <- unique(dt$gene_[!is.na(fv) & fv < fdr_thresh])
  } else {
    cat("Sem coluna de significancia/FDR em", f0, "\n")
    next
  }

  all_genes <- unique(dt$gene_)
  de_genes  <- de_genes[!is.na(de_genes)]
  cat("Genes no arquivo:", length(all_genes), "| DE:", length(de_genes), "\n")

  # Demais entra no namespace do catalogo para o Fisher (mesmo criterio do
  # cross_DEGs_STRs.py): se o arquivo de DEG usa simbolos, comparar por
  # gene_name; se usa ENSG, comparar por gene_id.
  deg_looks_ens <- (sum(grepl("^ENSG", all_genes)) / max(1, length(all_genes))) > 0.5
  str_genes <- if (deg_looks_ens) str_genes_ens else str_genes_sym

  if (length(de_genes) < min_de) {
    cat("Poucos genes DE (", length(de_genes), " < min-de ", min_de,
        ") -- pulando Fisher.\n", sep = "")
  } else {
    ff <- fisher_de_str(de_genes, all_genes, str_genes)
    ff <- cbind(data.frame(gse = gse), ff)
    res_fisher[[gse]] <- ff
    print(ff)
  }

  # overlap DE x genes com STR outlier (se houver arquivos de outlier no dir)
  outl_files <- list.files(d, pattern = "outlier", full.names = TRUE,
                           ignore.case = TRUE)
  if (length(outl_files)) {
    o <- tryCatch(fread(outl_files[1], header = TRUE, sep = "\t",
                        data.table = FALSE), error = function(e) NULL)
    if (!is.null(o)) {
      og_col <- detect_col(colnames(o), c("gene", "gene_name", "gene_symbol"))
      if (!is.na(og_col)) {
        og <- unique(as.character(o[[og_col]]))
        res_overlap[[gse]] <- data.frame(
          gse = gse, n_de = length(de_genes),
          n_de_with_str_outlier = length(intersect(de_genes, og)),
          de_with_str_outlier_genes = paste(intersect(de_genes, og),
                                            collapse = ";"))
      }
    }
  }
}

# ---------------------------------------------------------------------------
# 4. Saidas
# ---------------------------------------------------------------------------
fisher_tab <- bind_rows(res_fisher)
if (nrow(fisher_tab)) {
  fisher_tab$q_value <- p.adjust(fisher_tab$fisher_p, method = "BH")
  fwrite(fisher_tab, file.path(out_dir, "de_str_fisher_by_study.tsv"), sep = "\t")
  cat("\n=== Fisher DE x STR por estudo ===\n")
  print(fisher_tab)
} else {
  cat("Nenhum estudo com Fisher calculado.\n")
}

if (length(res_overlap)) {
  fwrite(bind_rows(res_overlap),
         file.path(out_dir, "de_outlier_str_overlap_by_study.tsv"), sep = "\t")
  cat("\n=== Overlap DE x STR-outlier por estudo ===\n")
  print(bind_rows(res_overlap))
}

# --- Enriquecimento de vias (hipergeometrico) nos genes DE (pool de estudos)
if (length(res_fisher)) {
  de_pool <- unique(unlist(lapply(deg_dirs, function(d) {
    f0 <- list.files(d, pattern = "(tsv|txt|csv)$", full.names = TRUE,
                     ignore.case = TRUE)[1]
    if (is.na(f0)) return(character(0))
    dt <- tryCatch(fread(f0, header = TRUE, sep = "\t", data.table = FALSE),
                   error = function(e) NULL)
    if (is.null(dt)) return(character(0))
    gene_col <- detect_col(colnames(dt), c("Gene", "gene_symbol", "gene_name"))
    sig_col  <- detect_col(colnames(dt), c("Significant", "significance"))
    fdr_col  <- detect_col(colnames(dt), c("FDR", "adj.P.Val", "P.Value"))
    if (is.na(gene_col)) return(character(0))
    g <- as.character(dt[[gene_col]])
    if (!is.na(sig_col)) {
      keep <- dt[[sig_col]] %in% c(1, TRUE, "TRUE", "yes", "Yes", "sig", "Sig", "Significant")
    } else if (!is.na(fdr_col)) {
      fv <- suppressWarnings(as.numeric(dt[[fdr_col]]))
      keep <- !is.na(fv) & fv < fdr_thresh
    } else return(character(0))
    g[keep & !is.na(g)]
  })))

  if (length(de_pool) >= min_de) {
    cat("\nPool de genes DE (todos estudos):", length(de_pool), "\n")
    # Converte para ENSG (espaco dos mapas de vias) via catalogo, se preciso.
    de_pool_ens <- if (any(grepl("^ENSG", de_pool)))
      unique(de_pool[grepl("^ENSG", de_pool)])
    else if (nrow(sym_map))
      unique(sym_map$gene_id[match(de_pool, sym_map$gene_name)]) else character(0)
    de_pool_ens <- unique(de_pool_ens[!is.na(de_pool_ens)])
    cat("Pool convertido p/ ENSG:", length(de_pool_ens), "\n")
    if (length(de_pool_ens) >= min_de) {
      kegg_map <- get_map("kegg", cache_dir)
      rt_map   <- get_map("reactome", cache_dir)
      kegg_e  <- enrich_hyper(de_pool_ens, str_genes_ens, kegg_map, "KEGG")
      react_e <- enrich_hyper(de_pool_ens, str_genes_ens, rt_map, "REACTOME")
      all_e   <- rbind(kegg_e, react_e)
      all_e$q_value <- p.adjust(all_e$p_value, method = "BH")
      fwrite(all_e, file.path(out_dir, "de_pathway_enrichment.tsv"), sep = "\t")
      sig_e <- all_e[!is.na(all_e$q_value) & all_e$q_value < 0.05, , drop = FALSE]
      cat("=== Vias enriquecidas em genes DE (q<0.05):", nrow(sig_e), "===\n")
      if (nrow(sig_e)) print(sig_e) else print(head(all_e, 10))
    } else {
      cat("Pool ENSG insuficiente (", length(de_pool_ens),
          ") -- pulando enriquecimento de vias.\n", sep = "")
    }
  }
}

cat("\n=== FIM 4_de_enrichment.R ===\n")