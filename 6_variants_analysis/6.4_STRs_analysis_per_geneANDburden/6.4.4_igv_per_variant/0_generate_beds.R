#!/usr/bin/env Rscript
# 0_generate_beds.R
# Gera BEDs + mapeamento BAM para IGV.js (6.4.4).
# Outputs sao escritos no dir de execucao (data/).
#
# Para cada STR com outlier em suggestive_strs_outliers.tsv:
#   - amostra COM variante  = outlier_samples
#   - amostra SEM variante  = outro sample_id do MESMO STRs_ID em
#       STRs_normalized_residuals.tsv, nao-outlier, com BAM indexado em
#       bam_dir, de menor |allele2_residuals| (mais proximo do esperado).

REPO_ROOT <- "/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis"

outlier_file <- file.path(REPO_ROOT,
  "6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.1_per_str_analysis/6.4.1.2_dbscan_subset_GWAS/results/suggestive_strs_outliers.tsv")
norm_file <- file.path(REPO_ROOT,
  "5_global_dbscan/norm_test/STRs_normalized_residuals.tsv")
bam_dir <- "/storage/users/tulio/Projeto_Luy_COVID/results/recal/"

out_with    <- "str_samples_with_variant.bed"
out_without <- "str_samples_without_variant.bed"
out_map     <- "str_samples_bams.tsv"

pad <- 0
set.seed(20260813)

stopifnot(file.exists(outlier_file))
stopifnot(file.exists(norm_file))
stopifnot(dir.exists(bam_dir))

split_samples <- function(s) {
  s <- trimws(as.character(s))
  unlist(strsplit(s, "[;,[:space:]]+"))
}

db <- read.delim(outlier_file, header = TRUE, stringsAsFactors = FALSE)
db$STRs_ID <- trimws(db$STRs_ID)
db <- db[nzchar(trimws(db$outlier_samples)), ]

cat("[0_generate_beds] STRs com outlier:", nrow(db), "\n")

variants <- data.frame(
  STRs_ID  = rep(db$STRs_ID, lengths(lapply(db$outlier_samples, split_samples))),
  variant  = unlist(lapply(db$outlier_samples, split_samples)),
  stringsAsFactors = FALSE)

id_parts <- strsplit(variants$STRs_ID, ":", fixed = TRUE)
variants$chr    <- vapply(id_parts, function(x) x[1], character(1))
variants$start1 <- as.integer(vapply(id_parts, function(x) x[2], character(1)))
variants$motif  <- vapply(id_parts, function(x) x[3], character(1))
variants$copy   <- as.integer(vapply(id_parts, function(x) x[4], character(1)))
variants$start0 <- variants$start1 - 1
variants$end    <- variants$start0 + nchar(variants$motif) * variants$copy
variants$gene   <- NA_character_

norm <- read.delim(norm_file, header = TRUE, stringsAsFactors = FALSE)
norm$STRs_ID   <- trimws(norm$STRs_ID)
norm$sample_id <- trimws(norm$sample_id)
norm$resid     <- suppressWarnings(as.numeric(norm$allele2_residuals))

gene_map <- norm[!is.na(norm$gene_name) & nzchar(norm$gene_name),
                 c("STRs_ID", "gene_name")]
gene_map <- gene_map[!duplicated(gene_map$STRs_ID), ]
variants$gene <- gene_map$gene_name[match(variants$STRs_ID, gene_map$STRs_ID)]
variants$gene[is.na(variants$gene)] <- "UNKNOWN"

find_bam <- function(bam_dir, sid) {
  cands <- unique(c(sid, paste0(sid, ".bam"), sub("\\.bam$", "", sid)))
  for (b in cands) {
    bfull   <- file.path(bam_dir, b)
    baifull <- sub("\\.bam$", ".bai", bfull)
    if (file.exists(bfull) && file.exists(baifull)) return(bfull)
  }
  NULL
}

clean_name <- function(sid) sub("\\.bam$", "", sid)

verify_tab <- data.frame(
  gene = character(), STRs_ID = character(),
  variant_sample = character(), variant_bam = character(), variant_indexed = logical(),
  control_sample = character(), control_bam = character(), control_indexed = logical(),
  stringsAsFactors = FALSE)

mapping <- data.frame(
  gene = character(), STRs_ID = character(), chr = character(),
  start0 = integer(), end = integer(),
  variant_sample = character(), variant_bam = character(),
  control_sample = character(), control_bam = character(),
  stringsAsFactors = FALSE)

with_rows    <- list()
without_rows <- list()

for (i in seq_len(nrow(variants))) {
  v_strs <- variants$STRs_ID[i]
  v_chr  <- variants$chr[i]
  v_s0   <- variants$start0[i]
  v_end  <- variants$end[i]
  v_gene <- variants$gene[i]

  v_sid_raw <- as.character(variants$variant[i])
  v_bam <- find_bam(bam_dir, v_sid_raw)
  v_indexed <- !is.null(v_bam)
  v_sample  <- clean_name(v_sid_raw)

  c_sample <- NA_character_
  c_bam    <- NULL
  sub <- norm[norm$STRs_ID == v_strs, ]
  if (nrow(sub) > 0) {
    cand <- sub[!is.na(sub$sample_id) &
                clean_name(sub$sample_id) != clean_name(v_sid_raw), ]
    cand <- cand[!is.na(cand$resid), ]
    if (nrow(cand) > 0) {
      cand$in_dir <- vapply(cand$sample_id, function(s)
        !is.null(find_bam(bam_dir, s)), logical(1))
      cand <- cand[cand$in_dir, ]
      if (nrow(cand) > 0) {
        cand <- cand[order(abs(cand$resid)), ]
        c_row  <- cand[1, ]
        c_bam  <- find_bam(bam_dir, c_row$sample_id)
        c_sample <- clean_name(c_row$sample_id)
      }
    }
  }
  c_indexed <- !is.null(c_bam)

  if (!v_indexed)
    warning("Sem BAM indexado p/ amostra-com-variante: ", v_sample, " (", v_strs, ")")
  if (!c_indexed)
    warning("Sem controle valido (BAM indexado) p/ locus: ", v_gene, " (", v_strs, ")")

  v_score <- 0
  c_score <- if (!is.na(c_sample)) {
    s <- norm[norm$STRs_ID == v_strs & clean_name(norm$sample_id) == c_sample, ]$resid
    s <- if (length(s) && !is.na(s[1])) round(abs(s[1]) * 100) else 0
  } else 0

  s0 <- max(0, v_s0 - pad)
  e0 <- v_end + pad

  if (v_indexed) {
    with_rows[[length(with_rows) + 1]] <- data.frame(
      chr = v_chr, start = s0, end = e0,
      name = paste0(v_strs, "_", v_sample),
      score = v_score, strand = ".", stringsAsFactors = FALSE)
  }
  if (c_indexed) {
    without_rows[[length(without_rows) + 1]] <- data.frame(
      chr = v_chr, start = s0, end = e0,
      name = paste0(v_strs, "_", c_sample),
      score = c_score, strand = ".", stringsAsFactors = FALSE)
  }

  verify_tab <- rbind(verify_tab, data.frame(
    gene = v_gene, STRs_ID = v_strs,
    variant_sample = if (v_indexed) v_sample else NA_character_,
    variant_bam = if (v_indexed) basename(v_bam) else NA_character_,
    variant_indexed = v_indexed,
    control_sample = if (c_indexed) c_sample else NA_character_,
    control_bam = if (c_indexed) basename(c_bam) else NA_character_,
    control_indexed = c_indexed,
    stringsAsFactors = FALSE))

  mapping <- rbind(mapping, data.frame(
    gene = v_gene, STRs_ID = v_strs, chr = v_chr,
    start0 = v_s0, end = v_end,
    variant_sample = if (v_indexed) v_sample else NA_character_,
    variant_bam = if (v_indexed) v_bam else NA_character_,
    control_sample = if (c_indexed) c_sample else NA_character_,
    control_bam = if (c_indexed) c_bam else NA_character_,
    stringsAsFactors = FALSE))
}

write_bed <- function(rows, out) {
  if (length(rows) == 0) {
    cat("Nenhuma linha valida p/", out, "- arquivo NAO escrito.\n")
    return(0)
  }
  df <- do.call(rbind, rows)
  df <- df[order(df$chr, df$start), ]
  write.table(df, out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Escrevi", nrow(df), "loci ->", out, "\n")
  nrow(df)
}

n_with    <- write_bed(with_rows, out_with)
n_without <- write_bed(without_rows, out_without)

write.table(mapping, out_map, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Escrevi", nrow(mapping), "loci ->", out_map, "\n")

cat("\n=== Resumo (verify_tab) ===\n")
print(verify_tab)
