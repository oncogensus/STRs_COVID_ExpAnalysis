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

strategy <- "gwas_burden"

if (strategy == "gwas_burden") {
  background_file <- file.path(REPO_ROOT,
    "6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.1_per_str_analysis/6.4.1.2_dbscan_subset_GWAS/results/suggestive_gene_strs.tsv")
  out_dir <- file.path(REPO_ROOT, "7_variants_analysis/burden_test/results_gwas_burden")
  remove_sample_outliers <- FALSE
} else if (strategy == "gwas") {
  outlier_file <- file.path(REPO_ROOT,
    "6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.1_per_str_analysis/6.4.1.2_dbscan_subset_GWAS/results/suggestive_strs_outliers.tsv")
  outlier_col <- "outlier_samples"
  out_dir <- file.path(REPO_ROOT, "7_variants_analysis/burden_test/results_gwas")
  remove_sample_outliers <- FALSE
} else {
  outlier_file <- file.path(REPO_ROOT, "5_dbscan/outliers_search/results_dbscan/outliers_per_str.tsv")
  outlier_col <- "outlier_samples"
  out_dir <- file.path(REPO_ROOT, "7_variants_analysis/burden_test/results")
  remove_sample_outliers <- TRUE
}

case_label       <- "case"
control_label    <- "control"
high_thresh_q    <- 0.95
min_strs_per_gene <- 1
skat_kernel      <- "linear.weighted"
selftest         <- FALSE

normalize_ids <- function(ids) {
  ids <- as.character(ids)
  ids <- toupper(ids)
  ids <- sub("(?i)[._-]?\\d*BAM.*$", "", ids)
  ids <- sub("(?<=\\d)-[0-9]$", "", ids, perl = TRUE)
  ids <- sub("-0+([0-9]+)", "-\\1", ids)
  trimws(ids)
}

load_inputs <- function() {
  norm <- fread(norm_file, header = TRUE, sep = "\t", data.table = FALSE)
  required <- c("STRs_ID", "sample_id", "sample_id_clean", "group",
                "gene_id", "gene_name", "region", "allele2_est")
  if (!all(required %in% colnames(norm)))
    stop("norm_file deve conter: ", paste(required, collapse = ", "))

  id_map <- norm %>% distinct(sample_id, sample_id_clean)
  str_meta <- norm %>% distinct(STRs_ID, gene_id, gene_name, region)

  if (strategy == "gwas_burden") {
    if (!file.exists(background_file)) stop("background_file nao encontrado: ", background_file)
    bg <- fread(background_file, header = TRUE, sep = "\t", data.table = FALSE)
    bg_strs <- unique(bg$strs_id)
    cat("[gwas_burden] STRs no background:", length(bg_strs), "\n")

    dos <- norm[norm$STRs_ID %in% bg_strs, c("STRs_ID", "sample_id_clean", "allele2_est")]
    dos <- dos[!is.na(dos$allele2_est), ]

    cat("[gwas_burden] STRs com dosagem:", length(unique(dos$STRs_ID)), "\n")
    cat("[gwas_burden] Amostras com dosagem:", length(unique(dos$sample_id_clean)), "\n")

    dos_dt <- as.data.table(dos)
    wide <- dcast(dos_dt, sample_id_clean ~ STRs_ID, value.var = "allele2_est",
                  fun.aggregate = function(x) x[1])
    wide_df <- as.data.frame(wide)
    rownames(wide_df) <- wide_df$sample_id_clean
    wide_df$sample_id_clean <- NULL
    M_dosage <- as.matrix(wide_df)
    M_dosage[is.na(M_dosage)] <- 0

  } else if (strategy == "gwas") {
    if (!file.exists(outlier_file)) stop("outlier_file nao encontrado: ", outlier_file)
    dbs <- fread(outlier_file, header = TRUE, sep = "\t", data.table = FALSE)
    if (!all(c("STRs_ID", outlier_col) %in% colnames(dbs)))
      stop("outlier_file deve conter STRs_ID e '", outlier_col, "'")

    outlier_long <- data.frame(STRs_ID = character(0), sample_id_clean = character(0))
    for (i in seq_len(nrow(dbs))) {
      raw <- as.character(dbs[[outlier_col]][i])
      samps <- unlist(strsplit(raw, "[;,[:space:]]+"))
      samps <- samps[nzchar(trimws(samps))]
      if (length(samps) == 0) next
      m <- id_map$sample_id_clean[match(samps, id_map$sample_id)]
      if (all(is.na(m))) m <- id_map$sample_id_clean[match(samps, id_map$sample_id_clean)]
      m <- m[!is.na(m)]
      if (length(m) == 0) next
      outlier_long <- rbind(outlier_long,
                            data.frame(STRs_ID = rep(dbs$STRs_ID[i], length(m)),
                                       sample_id_clean = m, stringsAsFactors = FALSE))
    }
    if (nrow(outlier_long) == 0) stop("Nenhum outlier mapeado de outlier_file")

    all_samples <- unique(covar_tmp$sample_id_clean)
    all_strs <- unique(outlier_long$STRs_ID)
    M_dosage <- matrix(0L, nrow = length(all_samples), ncol = length(all_strs),
                       dimnames = list(all_samples, all_strs))
    M_dosage[cbind(outlier_long$sample_id_clean, outlier_long$STRs_ID)] <- 1L

  } else {
    if (!file.exists(outlier_file)) stop("outlier_file nao encontrado: ", outlier_file)
    dbs <- fread(outlier_file, header = TRUE, sep = "\t", data.table = FALSE)
    if (!all(c("STRs_ID", outlier_col) %in% colnames(dbs)))
      stop("outlier_file deve conter STRs_ID e '", outlier_col, "'")

    outlier_long <- data.frame(STRs_ID = character(0), sample_id_clean = character(0))
    for (i in seq_len(nrow(dbs))) {
      raw <- as.character(dbs[[outlier_col]][i])
      samps <- unlist(strsplit(raw, "[;,[:space:]]+"))
      samps <- samps[nzchar(trimws(samps))]
      if (length(samps) == 0) next
      m <- id_map$sample_id_clean[match(samps, id_map$sample_id)]
      if (all(is.na(m))) m <- id_map$sample_id_clean[match(samps, id_map$sample_id_clean)]
      m <- m[!is.na(m)]
      if (length(m) == 0) next
      outlier_long <- rbind(outlier_long,
                            data.frame(STRs_ID = rep(dbs$STRs_ID[i], length(m)),
                                       sample_id_clean = m, stringsAsFactors = FALSE))
    }
    if (nrow(outlier_long) == 0) stop("Nenhum outlier mapeado de outlier_file")

    all_samples <- unique(norm$sample_id_clean)
    all_strs <- unique(outlier_long$STRs_ID)
    M_dosage <- matrix(0L, nrow = length(all_samples), ncol = length(all_strs),
                       dimnames = list(all_samples, all_strs))
    M_dosage[cbind(outlier_long$sample_id_clean, outlier_long$STRs_ID)] <- 1L
  }

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
  covar <- covar %>% filter(!is.na(group))

  common_samples <- intersect(rownames(M_dosage), covar$sample_id_clean)
  M_dosage <- M_dosage[common_samples, , drop = FALSE]

  list(covar = covar, M_dosage = M_dosage, str_meta = str_meta)
}

run_skat_per_gene <- function(M, str_meta, covar, min_strs_per_gene, skat_kernel) {
  common <- intersect(rownames(M), covar$sample_id_clean)
  Mc <- M[common, , drop = FALSE]
  covc <- covar %>% filter(sample_id_clean %in% common) %>%
    arrange(match(sample_id_clean, rownames(Mc)))
  rownames(covc) <- covc$sample_id_clean
  covc <- covc[rownames(Mc), ]

  obj <- SKAT_Null_Model(group ~ age + sex + EV1 + EV2 + EV3,
                         data = covc, out_type = "D")

  genes <- unique(str_meta$gene_id[!is.na(str_meta$gene_id)])
  res <- list()
  for (g in genes) {
    strs_g <- str_meta$STRs_ID[str_meta$gene_id == g &
                               str_meta$STRs_ID %in% colnames(Mc)]
    if (length(strs_g) == 0) {
      res[[g]] <- data.frame(gene = g, gene_name = NA,
                             n_variants = 0, str_ids = "",
                             p_value = NA, note = "no STRs in matrix")
      next
    }
    cols <- Mc[, strs_g, drop = FALSE]
    keep <- which(apply(cols, 2, function(x) var(x) > 0))
    if (length(keep) < min_strs_per_gene) {
      res[[g]] <- data.frame(gene = g, gene_name = NA,
                             n_variants = length(keep),
                             str_ids = paste(strs_g, collapse = ";"),
                             p_value = NA, note = "fewer than min variants with variance")
      next
    }
    Z <- as.matrix(cols[, keep, drop = FALSE])
    r <- tryCatch(SKAT(Z, obj, kernel = skat_kernel),
                  error = function(e) NULL)
    gname <- unique(str_meta$gene_name[str_meta$gene_id == g])
    res[[g]] <- data.frame(gene = g,
                           gene_name = ifelse(length(gname) > 0, gname[1], NA),
                           n_variants = ncol(Z),
                           str_ids = paste(strs_g, collapse = ";"),
                           p_value = ifelse(is.null(r), NA, r$p.value),
                           note = ifelse(is.null(r), "SKAT error", ""))
  }
  out <- bind_rows(res)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

write_hits <- function(df, pcol, qcol, base, out_dir) {
  unc <- df[!is.na(df[[pcol]]) & df[[pcol]] < 0.05, , drop = FALSE]
  cor <- df[!is.na(df[[qcol]]) & df[[qcol]] < 0.05, , drop = FALSE]
  fwrite(unc, file.path(out_dir, paste0(base, "_hits_uncorrected.tsv")), sep = "\t")
  fwrite(cor, file.path(out_dir, paste0(base, "_hits_corrected.tsv")), sep = "\t")
  list(unc = unc, cor = cor)
}

if (selftest) {
  set.seed(1)
  n <- 300
  covar <- data.frame(
    sample_id_clean = paste0("S", 1:n),
    age = rnorm(n, 60, 10),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    group = rbinom(n, 1, 0.5),
    EV1 = rnorm(n), EV2 = rnorm(n), EV3 = rnorm(n),
    stringsAsFactors = FALSE)
  genes <- rep(paste0("G", 1:20), each = 5)
  str_meta <- data.frame(STRs_ID = paste0("STR", 1:100),
                         gene_id = genes, gene_name = genes,
                         region = "x", stringsAsFactors = FALSE)
  M <- matrix(rnorm(n * 100), nrow = n, ncol = 100,
              dimnames = list(covar$sample_id_clean, paste0("STR", 1:100)))
  cat("=== SELFTEST [", strategy, "] ===\n")
  skat <- run_skat_per_gene(M, str_meta, covar, min_strs_per_gene, skat_kernel)
  print(head(skat))
} else {
  cat("=== Strategy:", strategy, "===\n")
  cat("Lendo inputs...\n")
  inp <- load_inputs()
  cat("Matrix:", nrow(inp$M_dosage), "amostras x", ncol(inp$M_dosage), "STRs com dosagem\n")
  cat("  [debug] range allele2_est:", range(inp$M_dosage), "\n")
  cat("  [debug] n nao-zero:", sum(inp$M_dosage != 0), "\n")
  cat("  [debug] genes unicos:", length(unique(inp$str_meta$gene_id)), "\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat("Rodando SKAT por gene...\n")
  skat <- run_skat_per_gene(inp$M_dosage, inp$str_meta, inp$covar,
                            min_strs_per_gene, skat_kernel)
  fwrite(skat, file.path(out_dir, "skat_per_gene.tsv"), sep = "\t")
  sk <- write_hits(skat, "p_value", "q_value", "skat_per_gene", out_dir)

  cat("Resultados em:", out_dir, "\n")
  cat("\n=== SKAT por gene (top 20) ===\n")
  print(head(skat, 20))

  cat("\n=== SKAT: hits SEM correcao (p<0.05):", nrow(sk$unc), "===\n")
  print(sk$unc)
  cat("\n=== SKAT: hits COM correcao BH (q<0.05):", nrow(sk$cor), "===\n")
  print(sk$cor)
}
