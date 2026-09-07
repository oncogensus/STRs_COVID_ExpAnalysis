# ============================================================================
# 00_common.R -- funcoes compartilhadas do modulo 6.4.7 (SKAT-O + vias)
# ----------------------------------------------------------------------------
# Carrega inputs padrao, monta matrix de dosagem, ajusta o modelo nulo uma vez
# e roda SKAT-O (burden + SKAT combinados com rho otimo) por conjunto de STRs.
# Usado por 1_skat_o_per_gene.R e 2_skat_o_per_pathway.R.
#
# Estrategias de background:
#   gwas_burden  -> suggestive_gene_strs.tsv  (STRs de genes GWAS-sugestivos)
#   rna_burden   -> rna_gene_strs.tsv         (STRs de genes DEGs)
#   full         -> todos os STRs com dosagem no norm_file
#
# As colunas sao identicas ao burden_gwas.R (6.4.3), mantendo compatibilidade.
# ============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  if (!requireNamespace("SKAT", quietly = TRUE)) {
    install.packages("SKAT", repos = "https://cloud.r-project.org")
  }
  library(SKAT)
})

REPO_ROOT <- "/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"

norm_file  <- file.path(REPO_ROOT, "5_global_dbscan/norm_test/STRs_normalized_residuals.tsv")
pca_file   <- file.path(REPO_ROOT, "4_ancestry/EthSEQ_Results_3D/Report.PCAcoord")
pheno_file <- file.path(REPO_ROOT, "samples/samples_infos.csv")

case_label    <- "case"
control_label <- "control"

background_by_strategy <- function(strategy) {
  switch(strategy,
    gwas_burden = file.path(REPO_ROOT,
      "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/suggestive_gene_strs.tsv"),
    rna_burden  = file.path(REPO_ROOT,
      "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.2_dbscan_subset_RNA/6.4.2.2_RNA_matrix/results/rna_gene_strs.tsv"),
    full        = NA_character_,
    stop("Estrategia desconhecida: ", strategy,
         " (use gwas_burden | rna_burden | full)"))
}

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

normalize_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- toupper(ids)
  ids <- sub("(?i)[._-]?\\d*BAM.*$", "", ids)
  ids <- sub("(?<=\\d)-[0-9]$", "", ids, perl = TRUE)
  ids <- sub("-0+([0-9]+)", "-\\1", ids)
  trimws(ids)
}

# ---------------------------------------------------------------------------
# Carrega fenotipo (samples_infos.csv) + PCs de ancestralidade (Report.PCAcoord)
# e devolve covar com: sample_id_clean, age, sex, group (1/0), EV1:3.
# ---------------------------------------------------------------------------
load_covar <- function() {
  pheno <- fread(pheno_file, header = TRUE, sep = ",", data.table = FALSE)
  pheno_col <- grep("^sample$", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(pheno_col)) stop("pheno_file deve ter coluna 'sample'")
  pheno$sample_id_clean <- normalize_ids(pheno[[pheno_col]])
  pheno <- pheno %>% select(sample_id_clean, age, sex, group) %>%
    mutate(age = as.numeric(age), sex = as.factor(sex))

  pca <- read.delim(pca_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pca_col <- grep("^sample", colnames(pca), ignore.case = TRUE, value = TRUE)[1]
  colnames(pca)[colnames(pca) == pca_col] <- "sample_id"
  pca$sample_id_clean <- normalize_ids(pca$sample_id)
  pca <- pca %>% select(sample_id_clean, EV1, EV2, EV3)

  covar <- inner_join(pheno, pca, by = "sample_id_clean")
  covar$group <- ifelse(covar$group %in% case_label, 1,
                 ifelse(covar$group %in% control_label, 0, NA))
  if (any(is.na(covar$group))) warning("Alguns grupos NA: ", sum(is.na(covar$group)))
  covar <- covar %>% filter(!is.na(group)) %>% as.data.frame()
  covar
}

# ---------------------------------------------------------------------------
# Carrega o norm_file e monta a matrix de dosagem (amostras x STRs) segundo a
# estrategia. Devolve list(covar, M_dosage, str_meta).
# ---------------------------------------------------------------------------
load_inputs <- function(strategy, background_opt = NULL) {
  norm <- fread(norm_file, header = TRUE, sep = "\t", data.table = FALSE)
  required <- c("STRs_ID", "sample_id", "sample_id_clean", "group",
                "gene_id", "gene_name", "region", "allele2_est")
  if (!all(required %in% colnames(norm)))
    stop("norm_file deve conter: ", paste(required, collapse = ", "))

  str_meta <- norm %>% distinct(STRs_ID, gene_id, gene_name, region)

  covar <- load_covar()

  if (strategy == "full") {
    dos <- norm[!is.na(norm$allele2_est),
                c("STRs_ID", "sample_id_clean", "allele2_est")]
  } else {
    bg_file <- if (!is.null(background_opt)) background_opt
               else background_by_strategy(strategy)
    if (!file.exists(bg_file)) stop("background nao encontrado: ", bg_file)
    bg <- fread(bg_file, header = TRUE, sep = "\t", data.table = FALSE)
    if (!"strs_id" %in% colnames(bg)) stop("background_file deve conter coluna 'strs_id'")
    bg_strs <- unique(bg$strs_id)
    cat(sprintf("[%s] STRs no background: %d\n", strategy, length(bg_strs)))
    dos <- norm[norm$STRs_ID %in% bg_strs & !is.na(norm$allele2_est),
                c("STRs_ID", "sample_id_clean", "allele2_est")]
  }

  cat(sprintf("[%s] STRs com dosagem: %d | amostras com dosagem: %d\n",
              strategy, length(unique(dos$STRs_ID)),
              length(unique(dos$sample_id_clean))))

  dos_dt <- as.data.table(dos)
  wide <- dcast(dos_dt, sample_id_clean ~ STRs_ID, value.var = "allele2_est",
                fun.aggregate = function(x) x[1])
  wide_df <- as.data.frame(wide)
  rownames(wide_df) <- wide_df$sample_id_clean
  wide_df$sample_id_clean <- NULL
  M <- as.matrix(wide_df)
  M[is.na(M)] <- 0

  common <- intersect(rownames(M), covar$sample_id_clean)
  M <- M[common, , drop = FALSE]

  list(covar = covar, M_dosage = M, str_meta = str_meta)
}

# ---------------------------------------------------------------------------
# Ajusta o modelo nulo (case/control + idade + sexo + PCs) uma unica vez.
# ---------------------------------------------------------------------------
build_null_model <- function(covar, M) {
  common <- intersect(rownames(M), covar$sample_id_clean)
  covc <- covar %>% filter(sample_id_clean %in% common) %>%
    arrange(match(sample_id_clean, rownames(M)))
  rownames(covc) <- covc$sample_id_clean
  covc <- covc[rownames(M), ]

  obj <- SKAT_Null_Model(group ~ age + sex + EV1 + EV2 + EV3,
                         data = covc, out_type = "D")
  obj
}

# ---------------------------------------------------------------------------
# Roda SKAT-O num conjunto de STRs (filtra as presentes na matrix com
# variancia > 0). Devolve data.frame de 1 linha com p/o, p/skat, p/burden,
# rho otimo, n_variants e note.
#
# Compatibilidade de API do pacote SKAT:
#   * SKAT <= 2.1: funcao SKAT.O(Z, obj, kernel) -> p.value/p.value.SKAT/
#     p.value.burden/rho (retorno direto).
#   * SKAT >= 2.2: SKAT.O removida; usar SKAT(Z, obj, kernel,
#     r.corr = seq(0,1,by=0.1)) que retorna o p combinado em $p.value, os
#     p por rho em $param$p.val.each e o rho estimado em $param$rho_est.
#     Os betas rho=0 (SKAT) e rho=1 (burden) sao as 1a e ultima entradas.
# ---------------------------------------------------------------------------
run_skat_o <- function(str_ids, M, obj, min_var = 1, kernel = "linear.weighted") {
  str_ids <- as.character(str_ids)
  present <- str_ids[str_ids %in% colnames(M)]
  if (length(present) == 0) {
    return(data.frame(p_value_skat_o = NA_real_, p_value_skat = NA_real_,
                      p_value_burden = NA_real_, rho = NA_real_,
                      n_variants = 0L, note = "no STRs in matrix"))
  }
  Z0 <- as.matrix(M[, present, drop = FALSE])
  keep <- which(apply(Z0, 2, function(x) var(x) > 0))
  if (length(keep) < min_var) {
    return(data.frame(p_value_skat_o = NA_real_, p_value_skat = NA_real_,
                      p_value_burden = NA_real_, rho = NA_real_,
                      n_variants = length(keep),
                      note = "fewer than min variants with variance"))
  }
  Z <- as.matrix(Z0[, keep, drop = FALSE])

  skat_o_old <- function(Zm) {
    r <- tryCatch(SKAT.O(Zm, obj, kernel = kernel), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    list(p_value_skat_o = r$p.value, p_value_skat = r$p.value.SKAT,
         p_value_burden = r$p.value.burden, rho = r$rho, note = "")
  }
  skat_o_new <- function(Zm) {
    rho_grid <- seq(0, 1, by = 0.1)
    r <- tryCatch(SKAT(Zm, obj, kernel = kernel, r.corr = rho_grid),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    p_each <- r$param$p.val.each
    if (is.null(p_each) || !length(p_each)) return(NULL)
    n <- length(p_each)
    rho_est <- if (!is.null(r$param$rho_est)) r$param$rho_est
               else rho_grid[which.min(p_each)]
    list(p_value_skat_o = r$p.value,
         p_value_skat   = p_each[1],
         p_value_burden = p_each[n],
         rho            = rho_est,
         note           = "api: SKAT r.corr-grid (SKAT.O ausente)")
  }

  r <- skat_o_old(Z)
  if (is.null(r)) r <- skat_o_new(Z)
  if (is.null(r)) {
    single <- tryCatch(SKAT(Z, obj, kernel = kernel), error = function(e) NULL)
    if (is.null(single)) {
      return(data.frame(p_value_skat_o = NA_real_, p_value_skat = NA_real_,
                        p_value_burden = NA_real_, rho = NA_real_,
                        n_variants = ncol(Z), note = "SKAT.O/SKAT error"))
    }
    return(data.frame(p_value_skat_o = single$p.value,
                      p_value_skat   = single$p.value,
                      p_value_burden = single$p.value,
                      rho            = NA_real_,
                      n_variants     = ncol(Z),
                      note           = "fallback: SKAT unico"))
  }
  data.frame(p_value_skat_o = r$p_value_skat_o,
             p_value_skat   = r$p_value_skat,
             p_value_burden = r$p_value_burden,
             rho            = r$rho,
             n_variants     = ncol(Z),
             note           = r$note)
}

# ---------------------------------------------------------------------------
# Grava tabela completa + hits uncorrected/corrected (BH).
# ---------------------------------------------------------------------------
write_hits <- function(df, pcol, qcol, base, out_dir) {
  fwrite(df, file.path(out_dir, paste0(base, ".tsv")), sep = "\t")
  unc <- df[!is.na(df[[pcol]]) & df[[pcol]] < 0.05, , drop = FALSE]
  cor <- df[!is.na(df[[qcol]]) & df[[qcol]] < 0.05, , drop = FALSE]
  fwrite(unc, file.path(out_dir, paste0(base, "_hits_uncorrected.tsv")), sep = "\t")
  fwrite(cor, file.path(out_dir, paste0(base, "_hits_corrected.tsv")), sep = "\t")
  list(unc = unc, cor = cor)
}

# ---------------------------------------------------------------------------
# Mapas de vias (msigdbr C2 CP:KEGG / CP:REACTOME) com cache em TSV.
# ---------------------------------------------------------------------------
load_pathway_map_msigdbr <- function(subcat) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    cat("Aviso: msigdbr ausente (mapa vazio).\n")
    return(data.frame(pathway = character(0), pathway_name = character(0),
                      gene_id = character(0), stringsAsFactors = FALSE))
  }
  gs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                         subcategory = subcat)
  data.frame(pathway = as.character(gs$gs_name),
             pathway_name = as.character(gs$gs_name),
             gene_id = as.character(gs$ensembl_gene),
             stringsAsFactors = FALSE)
}

get_map <- function(tag, cache_dir) {
  cache_file <- file.path(cache_dir, paste0(tag, "_pathway_map.tsv"))
  if (file.exists(cache_file)) {
    cat("USANDO mapa cacheado:", cache_file, "\n")
    return(fread(cache_file, header = TRUE, sep = "\t", data.table = FALSE)[,
           c("pathway", "pathway_name", "gene_id")])
  }
  m <- load_pathway_map_msigdbr(if (tag == "kegg") "CP:KEGG" else "CP:REACTOME")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(m, cache_file, sep = "\t")
  m
}

# ---------------------------------------------------------------------------
# Teste hipergeometrico (phyper) de enriquecimento em vias no espaco Ensembl.
# retorna tibble ordenado por p_value, com BH.
# ---------------------------------------------------------------------------
enrich_hyper <- function(genes_i, universe_i, map, base) {
  map <- map[map$gene_id %in% c(genes_i, universe_i), , drop = FALSE]
  if (!nrow(map)) return(NULL)
  interes  <- intersect(unique(genes_i), unique(map$gene_id))
  universe_i <- if (length(universe_i)) unique(universe_i) else unique(map$gene_id)
  univ_eff <- intersect(unique(universe_i), unique(map$gene_id))
  N <- length(unique(univ_eff))
  if (N == 0) return(NULL)
  k_i <- length(intersect(interes, univ_eff))
  rows <- split(map$gene_id, map$pathway)
  res <- lapply(unique(map$pathway), function(pw) {
    pw_genes <- intersect(unique(rows[[pw]]), univ_eff)
    q <- length(intersect(pw_genes, interes))
    m <- length(pw_genes)
    pv <- if (q == 0) 1 else
      phyper(q - 1, m, N - m, k_i, lower.tail = FALSE)
    data.frame(base = base,
               pathway = pw,
               pathway_name = unique(map$pathway_name[map$pathway == pw])[1],
               n_genes_in_pathway = m,
               n_overlap = q,
               genes_overlap = paste(intersect(pw_genes, interes), collapse = ";"),
               p_value = pv, q_value = NA_real_)
  })
  out <- bind_rows(res)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}