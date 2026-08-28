#!/usr/bin/env Rscript
# run_dbscan_subset.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Re-roda o DBSCAN (mesmos parametros de 5_global_dbscan/outliers_search/5.2_dbscan_str.r)
#   sobre o SUBSET de STRs localizados nos genes sugestivos do COVID-19 HG r7
#   (p < 1e-5). Para cada STR do subset, detecta amostras "outlier" nos residuos
#   normalizados de alelo2 (allele2_residuals).
#
# DEFINICAO DE OUTLIER (AD_STR-style)
#   - DBSCAN 1D nos residuos de cada STR: pontos em clusters = "normais", pontos
#     em ruido (cluster 0) = candidatos a outlier.
#   - cutoff = maior residual ENTRE os pontos nao-ruido (clusters != 0), com piso
#     de 2 (para o residuo ser considerado biologicamente relevante).
#   - Outlier = amostra com residual > cutoff.
#   - Se o STR nao tem ruido OU tem 1 unico cluster (todos agrupados), cutoff = Inf
#     -> nenhum outlier (nao ha separacao).
#   - STRs com < 3 residuos validos sao pulados (DBSCAN instavel).
#
# PARAMETROS DBSCAN (heuristicos, dependem de n_valid)
#   minPts = max(2, ceiling(log2(2 * n_valid)))   # cresce com log2 do n de amostras
#   eps     = IQR(residuos) = quantile(0.95)-quantile(0.05), piso 1e-6
#
# DEPENDENCIAS (env micromamba `dbscan-r`)
#   micromamba install -n dbscan-r -c conda-forge -y r-data.table r-dbscan
#
# USO
#   Rscript run_dbscan_subset.R <input_residuals.tsv> <output.tsv> [verbose]
#     verbose: "true"/"1" ativa debug por STR (opcional)
# ---------------------------------------------------------------------------
for (pkg in c("data.table", "dbscan")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Pacote '", pkg, "' ausente no env `dbscan-r`. Instale com:\n",
         "  micromamba install -n dbscan-r -c conda-forge -y r-", pkg, "\n")
  }
}
library(data.table)
library(dbscan)

# --- argumentos ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript run_dbscan_subset.R <input> <output> [verbose]")
input_file  <- args[1]
output_file <- args[2]
verbose <- length(args) >= 3 && tolower(args[3]) %in% c("true", "1", "yes", "verbose")
cat(sprintf("[DEBUG] input=%s output=%s verbose=%s\n", input_file, output_file, verbose))

# --- leitura ---
minpts <- 2
dat <- fread(input_file, header = TRUE, sep = "\t", data.table = FALSE)
cat(sprintf("[DEBUG] linhas lidas=%d  colunas=%d\n", nrow(dat), ncol(dat)))

# 'group' e reservado (vem do arquivo de residuos) mas nao usado neste script.
required_cols <- c("STRs_ID", "allele2_residuals", "sample_id", "group")
if (!all(required_cols %in% colnames(dat)))
  stop("Input deve conter: ", paste(required_cols, collapse = ", "))

# --- lista de STRs unicos ---
str_list <- unique(dat$STRs_ID)
cat(sprintf("Total STRs (subset): %d\n", length(str_list)))

# tabela de resultados (1 linha por STR)
results <- data.frame(
  STRs_ID          = str_list,
  n_samples        = NA,   # total de amostras (com/sem residuo)
  n_samples_valid  = NA,   # amostras com residuo nao-NA
  n_outliers       = NA,
  outlier_samples  = NA,   # amostras outlier (separadas por ;)
  outlier_residuals = NA,  # residuos das amostras outlier (separados por ;)
  n_clusters       = NA,   # n de clusters nao-ruido
  noise_ratio      = NA,   # proporcao de ruido
  eps              = NA_real_,  # raio DBSCAN (IQR dos residuos, piso 1e-6)
  minPts           = NA_integer_,# minPts DBSCAN (cresce com log2(n_valid), piso 2)
  cutoff           = NA_real_,  # limiar de residual p/ outlier (maior residual nao-ruido, piso 2)
  max_residual     = NA_real_,  # maior |residuo| do STR (contexto dos residuos)
  mean_residual    = NA_real_,  # residuo medio do STR
  stringsAsFactors = FALSE
)

# contadores para o resumo final
n_skipped_lowN  <- 0   # STRs com < 3 residuos validos
n_with_outliers <- 0   # STRs com >= 1 outlier
total_outliers  <- 0   # total de amostras-outlier

for (i in seq_along(str_list)) {
  if (i %% 100 == 0) cat(sprintf("Processing STR %d of %d\n", i, length(str_list)))

  str_id  <- str_list[i]
  dat_str <- subset(dat, STRs_ID == str_id)
  n_total <- nrow(dat_str)
  results[i, "n_samples"] <- n_total

  # so usa amostras com residuo valido
  dat_str <- dat_str[!is.na(dat_str$allele2_residuals), ]
  n_valid <- nrow(dat_str)
  results[i, "n_samples_valid"] <- n_valid

  # DBSCAN instavel com poucas amostras -> pula (sem outliers)
  if (n_valid < 3) {
    n_skipped_lowN <- n_skipped_lowN + 1
    if (verbose) cat(sprintf("  [%s] pulado: n_valid=%d (<3)\n", str_id, n_valid))
    next
  }
  residuals <- dat_str$allele2_residuals
  sample_ids <- dat_str$sample_id

  # minPts cresce com log2 do n de amostras (piso 2)
  minPts <- ceiling(log2(minpts * n_valid))
  minPts <- max(2, minPts)
  # eps = IQR dos residuos (piso 1e-6 para evitar eps zero)
  eps_value <- as.numeric(diff(quantile(residuals, c(0.05, 0.95), na.rm = TRUE)))
  eps_value <- max(eps_value, 1e-6)
  if (verbose) cat(sprintf("  [%s] n_valid=%d minPts=%d eps=%.4g\n",
                           str_id, n_valid, minPts, eps_value))

  # DBSCAN 1D (cada STR e 1 dimensao de residuos)
  X <- matrix(residuals, ncol = 1)
  scan <- dbscan(X, eps = eps_value, minPts = minPts)
  clusters <- scan$cluster
  unique_clusters <- unique(clusters)

  # cutoff de outlier:
  #  - se nao ha ruido (todos agrupados) OU so 1 cluster -> sem separacao -> Inf
  #  - caso contrario, maior residual dos pontos nao-ruido, com piso 2
  if (length(unique_clusters) == 1 || sum(clusters == 0) == 0) {
    cutoff <- Inf
  } else {
    non_noise_residuals <- residuals[clusters != 0]
    cutoff <- max(non_noise_residuals, na.rm = TRUE)
    cutoff <- ifelse(cutoff < 2, 2, cutoff)
  }
  if (verbose) cat(sprintf("  [%s] n_clusters=%d noise_ratio=%.3f cutoff=%.3g\n",
                           str_id, length(unique_clusters[unique_clusters != 0]),
                           sum(clusters == 0) / n_valid, cutoff))

  # amostras com residual estritamente acima do cutoff sao outliers
  outlier_idx <- which(residuals > cutoff)
  n_outliers <- length(outlier_idx)
  results[i, "n_outliers"] <- n_outliers
  results[i, "outlier_samples"]   <- paste(sample_ids[outlier_idx], collapse = ";")
  results[i, "outlier_residuals"] <- paste(residuals[outlier_idx], collapse = ";")
  results[i, "n_clusters"] <- length(unique_clusters[unique_clusters != 0])
  results[i, "noise_ratio"] <- sum(clusters == 0) / n_valid
  results[i, "eps"] <- eps_value
  results[i, "minPts"] <- minPts
  # cutoff Inf indica "sem separacao" -> NA para nao poluir o TSV com Inf
  results[i, "cutoff"] <- ifelse(is.finite(cutoff), cutoff, NA_real_)
  results[i, "max_residual"] <- max(residuals, na.rm = TRUE)
  results[i, "mean_residual"] <- mean(residuals, na.rm = TRUE)

  # acumula contadores do resumo
  if (n_outliers > 0) {
    n_with_outliers <- n_with_outliers + 1
    total_outliers  <- total_outliers + n_outliers
  }
}

# --- escrita + resumo ---
fwrite(results, output_file, sep = "\t", row.names = FALSE)
cat(sprintf("Results saved to: %s\n", output_file))
cat(sprintf("[RESUMO] STRs=%d | pulados(n_valid<3)=%d | com outlier=%d | amostras-outlier=%d\n",
            length(str_list), n_skipped_lowN, n_with_outliers, total_outliers))
