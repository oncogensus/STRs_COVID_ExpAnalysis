suppressMessages({
  library(data.table)
  library(dplyr)
  if (!requireNamespace("SKAT", quietly = TRUE)) {
    install.packages("SKAT", repos = "https://cloud.r-project.org")
  }
  library(SKAT)
})

REPO_ROOT  <- "/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"

norm_file  <- file.path(REPO_ROOT, "5_global_dbscan/norm_test/STRs_normalized_residuals.tsv")
pca_file   <- file.path(REPO_ROOT, "4_ancestry/EthSEQ_Results_3D/Report.PCAcoord")
pheno_file <- file.path(REPO_ROOT, "samples/samples_infos.csv")

strategy <- "gwas"

if (strategy == "gwas") {
  outlier_file     <- file.path(REPO_ROOT, "6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.1_per_str_analysis/6.4.1.2_dbscan_subset_GWAS/results/suggestive_strs_outliers.tsv")
  outlier_col      <- "outlier_samples"
  out_dir          <- file.path(REPO_ROOT, "7_variants_analysis/burden_test/results_gwas")
  remove_sample_outliers <- FALSE
} else {
  outlier_file     <- file.path(REPO_ROOT, "5_dbscan/outliers_search/results_dbscan/outliers_per_str.tsv")
  outlier_col      <- "outlier_samples"
  out_dir          <- file.path(REPO_ROOT, "7_variants_analysis/burden_test/results")
  remove_sample_outliers <- TRUE
}

case_label           <- "case"
control_label        <- "control"
high_thresh_q        <- 0.95
min_strs_per_gene    <- 2
skat_kernel          <- "linear.weighted"
run_gene_burden      <- TRUE
selftest             <- FALSE

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
                "gene_id", "gene_name", "region")
  if (!all(required %in% colnames(norm)))
    stop("norm_file deve conter: ", paste(required, collapse = ", "))

  id_map <- norm %>% distinct(sample_id, sample_id_clean)

  if (!file.exists(outlier_file)) stop("outlier_file nao encontrado: ", outlier_file,
    "\nAjuste 'outlier_file' ou 'strategy'.")
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

  str_meta <- norm %>% distinct(STRs_ID, gene_id, gene_name, region)

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
  if (any(is.na(covar$group))) warning("Alguns grupos NA (nem case nem control): ",
                                        sum(is.na(covar$group)))
  covar <- covar %>% filter(!is.na(group))

  list(covar = covar, outlier_long = outlier_long, str_meta = str_meta)
}

build_matrix <- function(covar, outlier_long) {
  all_samples <- covar$sample_id_clean
  all_strs <- unique(outlier_long$STRs_ID)
  M <- matrix(0, nrow = length(all_samples), ncol = length(all_strs),
              dimnames = list(all_samples, all_strs))
  M[cbind(outlier_long$sample_id_clean, outlier_long$STRs_ID)] <- 1
  M
}

run_burden_global <- function(M, covar, high_thresh_q) {
  count <- rowSums(M)
  cat("  [burden_global] range(count)=", range(count),
      " n_zero=", sum(count == 0), " n_nonzero=", sum(count > 0), "\n")
  cat("  [burden_global] table(group)=", table(covar$group[match(rownames(M), covar$sample_id_clean)]), "\n")

  df <- data.frame(
    sample_id_clean = rownames(M),
    count = count,
    group = covar$group[match(rownames(M), covar$sample_id_clean)],
    age = covar$age[match(rownames(M), covar$sample_id_clean)],
    sex = covar$sex[match(rownames(M), covar$sample_id_clean)],
    EV1 = covar$EV1[match(rownames(M), covar$sample_id_clean)],
    EV2 = covar$EV2[match(rownames(M), covar$sample_id_clean)],
    EV3 = covar$EV3[match(rownames(M), covar$sample_id_clean)],
    stringsAsFactors = FALSE
  )

  empty_row <- data.frame(
    test = c("logistic_OR_per_expansion", "mann_whitney", "threshold_OR", "threshold_fisher"),
    estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, stringsAsFactors = FALSE)

  if (var(count) == 0 || length(unique(df$group)) < 2) {
    cat("  [burden_global] count sem variancia ou grupo unico: retornando NA\n")
    return(empty_row)
  }

  fit <- tryCatch(
    glm(group ~ count + age + sex + EV1 + EV2 + EV3,
        data = df, family = binomial()),
    error = function(e) { cat("  [burden_global] glm erro:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(empty_row)
  coefs <- tryCatch(summary(fit)$coefficients["count", , drop = FALSE],
                    error = function(e) NULL)
  if (is.null(coefs) || !"count" %in% rownames(coefs)) {
    cat("  [burden_global] 'count' nao estimavel (modelo colapsou): retornando NA\n")
    return(empty_row)
  }
  or <- exp(coefs[, "Estimate"])
  or_lo <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  or_hi <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])

  w <- wilcox.test(count ~ group, data = df)

  hi <- quantile(df$count, high_thresh_q)
  tbl <- table(df$count >= hi, df$group)
  fisher_p <- if (all(dim(tbl) == c(2, 2))) fisher.test(tbl)$p.value else NA
  or_hi_thr <- if (all(dim(tbl) == c(2, 2))) {
    a <- tbl[2, 2]; b <- tbl[2, 1]; c <- tbl[1, 2]; d <- tbl[1, 1]
    (a / (a + b)) / (c / (c + d))
  } else NA

  data.frame(
    test = c("logistic_OR_per_expansion", "mann_whitney", "threshold_OR", "threshold_fisher"),
    estimate = c(or, NA, or_hi_thr, NA),
    ci_low = c(or_lo, NA, NA, NA),
    ci_high = c(or_hi, NA, NA, NA),
    p_value = c(coefs[, "Pr(>|z|)"], w$p.value, NA, fisher_p),
    stringsAsFactors = FALSE
  )
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
    cols <- Mc[, strs_g, drop = FALSE]
    keep <- which(apply(cols, 2, function(x) var(x) > 0))
    if (length(keep) < min_strs_per_gene) {
      res[[g]] <- data.frame(gene = g, gene_name = NA,
                             n_variants = length(keep),
                             p_value = NA, note = "fewer than min variants")
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
  out
}

run_gene_burden_test <- function(M, str_meta, covar) {
  common <- intersect(rownames(M), covar$sample_id_clean)
  Mc <- M[common, , drop = FALSE]
  covc <- covar %>% filter(sample_id_clean %in% common) %>%
    arrange(match(sample_id_clean, rownames(Mc)))
  rownames(covc) <- covc$sample_id_clean
  covc <- covc[rownames(Mc), ]

  genes <- unique(str_meta$gene_id[!is.na(str_meta$gene_id)])
  res <- list()
  for (g in genes) {
    strs_g <- str_meta$STRs_ID[str_meta$gene_id == g &
                               str_meta$STRs_ID %in% colnames(Mc)]
    if (length(strs_g) == 0) next
    gcount <- rowSums(Mc[, strs_g, drop = FALSE])
    if (var(gcount) == 0) next
    df <- data.frame(
      count = gcount, group = covc$group, age = covc$age,
      sex = covc$sex, EV1 = covc$EV1, EV2 = covc$EV2, EV3 = covc$EV3,
      stringsAsFactors = FALSE)
    fit <- glm(group ~ count + age + sex + EV1 + EV2 + EV3,
               data = df, family = binomial())
    cf <- summary(fit)$coefficients["count", ]
    gname <- unique(str_meta$gene_name[str_meta$gene_id == g])
    res[[g]] <- data.frame(gene = g,
                           gene_name = ifelse(length(gname) > 0, gname[1], NA),
                           n_variants = length(strs_g),
                           str_ids = paste(strs_g, collapse = ";"),
                           p_value = cf["Pr(>|z|)"],
                           OR_per_expansion = exp(cf["Estimate"]),
                           stringsAsFactors = FALSE)
  }
  out <- bind_rows(res)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

run_str_burden_test <- function(M, covar) {
  common <- intersect(rownames(M), covar$sample_id_clean)
  Mc <- M[common, , drop = FALSE]
  covc <- covar %>% filter(sample_id_clean %in% common) %>%
    arrange(match(sample_id_clean, rownames(Mc)))
  rownames(covc) <- covc$sample_id_clean
  covc <- covc[rownames(Mc), ]

  res <- list()
  for (s in colnames(Mc)) {
    if (var(Mc[, s]) == 0) next
    df <- data.frame(
      str = Mc[, s], group = covc$group, age = covc$age,
      sex = covc$sex, EV1 = covc$EV1, EV2 = covc$EV2, EV3 = covc$EV3,
      stringsAsFactors = FALSE)
    fit <- tryCatch(glm(group ~ str + age + sex + EV1 + EV2 + EV3,
                        data = df, family = binomial()),
                    error = function(e) NULL)
    if (is.null(fit)) next
    cf <- tryCatch(summary(fit)$coefficients["str", ], error = function(e) NULL)
    if (is.null(cf)) next
    res[[s]] <- data.frame(STRs_ID = s,
                           n_outliers = sum(Mc[, s] == 1),
                           p_value = cf["Pr(>|z|)"],
                           OR_per_expansion = exp(cf["Estimate"]),
                           stringsAsFactors = FALSE)
  }
  out <- bind_rows(res)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

run_pipeline <- function(covar, outlier_long, str_meta) {
  if (remove_sample_outliers && nrow(outlier_long) > 0) {
    tmp <- build_matrix(covar, outlier_long)
    cnt <- rowSums(tmp)
    thr <- mean(cnt) + 2 * sd(cnt)
    bad <- rownames(tmp)[cnt > thr]
    if (length(bad) > 0) {
      cat("Removendo", length(bad), "amostras outlier (>2 DP de expansoes):",
          paste(bad, collapse = ", "), "\n")
      covar <- covar %>% filter(!sample_id_clean %in% bad)
      outlier_long <- outlier_long %>% filter(!sample_id_clean %in% bad)
    }
  }
  M <- build_matrix(covar, outlier_long)
  cat("Matriz [", strategy, "]:", nrow(M), "amostras x", ncol(M), "STRs com expansao\n")
  cat("  [debug] outlier_long nrow =", nrow(outlier_long), "\n")
  cat("  [debug] n de amostras unicas em outlier_long =",
      length(unique(outlier_long$sample_id_clean)), "\n")
  cat("  [debug] contagens por STR (colSums M):\n")
  print(colSums(M))
  cat("  [debug] head(outlier_long):\n")
  print(head(outlier_long, 12))

  burden <- run_burden_global(M, covar, high_thresh_q)
  skat <- run_skat_per_gene(M, str_meta, covar, min_strs_per_gene, skat_kernel)
  out <- list(burden = burden, skat = skat)
  if (run_gene_burden) out$gene_burden <- run_gene_burden_test(M, str_meta, covar)
  out$str_burden <- run_str_burden_test(M, covar)
  out
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
  outs <- sample(1:100, 40)
  outlier_long <- data.frame(
    STRs_ID = rep(paste0("STR", outs), each = 3),
    sample_id_clean = sample(covar$sample_id_clean, 120, TRUE),
    stringsAsFactors = FALSE) %>% distinct()
  cat("=== SELFTEST [", strategy, "] ===\n")
  res <- run_pipeline(covar, outlier_long, str_meta)
  print(res$burden)
  print(head(res$skat[order(res$skat$p_value), ]))
  if (!is.null(res$gene_burden)) print(head(res$gene_burden[order(res$gene_burden$p_value), ]))
} else {
  inp <- load_inputs()
  res <- run_pipeline(inp$covar, inp$outlier_long, inp$str_meta)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(res$burden, file.path(out_dir, "burden_global.tsv"), sep = "\t")
  fwrite(res$skat, file.path(out_dir, "skat_per_gene.tsv"), sep = "\t")
  if (!is.null(res$gene_burden))
    fwrite(res$gene_burden, file.path(out_dir, "gene_burden.tsv"), sep = "\t")
  if (!is.null(res$str_burden))
    fwrite(res$str_burden, file.path(out_dir, "str_burden.tsv"), sep = "\t")

  sk <- write_hits(res$skat, "p_value", "q_value", "skat_per_gene", out_dir)
  if (!is.null(res$gene_burden))
    gb <- write_hits(res$gene_burden, "p_value", "q_value", "gene_burden", out_dir)
  sb <- write_hits(res$str_burden, "p_value", "q_value", "str_burden", out_dir)

  cat("Resultados em:", out_dir, "\n")
  cat("Burden global:\n"); print(res$burden)

  cat("\n=== SKAT por gene: hits SEM correção (p<0.05):", nrow(sk$unc), "===\n")
  print(sk$unc)
  cat("\n=== SKAT por gene: hits COM correção BH (q<0.05):", nrow(sk$cor), "===\n")
  print(sk$cor)

  if (!is.null(res$gene_burden)) {
    cat("\n=== Gene burden: hits SEM correção (p<0.05):", nrow(gb$unc), "===\n")
    print(gb$unc)
    cat("\n=== Gene burden: hits COM correção BH (q<0.05):", nrow(gb$cor), "===\n")
    print(gb$cor)
  }

  cat("\n=== STR burden (por STR): hits SEM correção (p<0.05):", nrow(sb$unc), "===\n")
  print(sb$unc)
  cat("\n=== STR burden (por STR): hits COM correção BH (q<0.05):", nrow(sb$cor), "===\n")
  print(sb$cor)
}
