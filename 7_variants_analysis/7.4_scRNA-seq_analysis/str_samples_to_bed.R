#!/usr/bin/env Rscript
# str_samples_to_bed.R
# Gera BEDs das amostras COM e SEM a variante, reutilizando a logica da
# celula 8 do 7.4.3.2_trackviewer_STR_viz.ipynb (alinhamento BAM variante vs controle).
# Standalone: NAO depende de trackViewer / GenomicRanges / Gviz.
# Rodar a partir de 7.4_scRNA-seq_analysis/ (onde STR_variants_UCSC_track.bed esta).
#
# Para cada um dos 8 loci de STR:
#   - amostra COM variante  = coluna `outlier_samples` do CSV unificado
#   - amostra SEM variante  = outro sample_id do MESMO STRs_ID em
#       STRs_normalized_residuals.tsv, nao-outlier, com BAM+`.bai` indexado em
#       bam_dir, de menor |allele2_residuals| (mais proximo do esperado).
# name do BED = paste0(STRs_ID, "_", sample)  (formato variante_nomeamostra).

# ---- inputs (edite conforme o ambiente) ----
variants_file <- '../7.3_STRs_filter/results/STR_vs_scRNA_overlap_unified.csv'
bed_file      <- 'STR_variants_UCSC_track.bed'
norm_file     <- '/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis/5_dbscan/norm_test/STRs_normalized_residuals.tsv'
bam_dir       <- '/storage/users/tulio/Projeto_Luy_COVID/results/recal/'

out_with      <- 'str_samples_with_variant.bed'
out_without   <- 'str_samples_without_variant.bed'

pad           <- 0      # 0 = janela exata do STR (do UCSC BED); use 20000 p/ +-20kb
set.seed(20260813)

stopifnot(file.exists(variants_file))
stopifnot(file.exists(norm_file))
stopifnot(dir.exists(bam_dir))

# ---- leitura do CSV de variantes ----
if (requireNamespace('readr', quietly = TRUE)) {
  ov <- readr::read_csv(variants_file, show_col_types = FALSE)
} else {
  ov <- read.csv(variants_file, stringsAsFactors = FALSE)
}
variants <- ov[!duplicated(ov$STRs_ID), ]

# STRs_ID codifica chr:pos:motif:copy (pos = 1-based start)
id_parts <- strsplit(variants$STRs_ID, ':', fixed = TRUE)
variants$chr    <- vapply(id_parts, function(x) x[1], character(1))
variants$start1 <- as.integer(vapply(id_parts, function(x) x[2], character(1)))
variants$motif  <- vapply(id_parts, function(x) x[3], character(1))
variants$copy   <- as.integer(vapply(id_parts, function(x) x[4], character(1)))
variants$start0 <- variants$start1 - 1
variants$gene   <- variants$gene_name

# end: do BED do UCSC se existir, senao estima por motif*copy
if (file.exists(bed_file)) {
  bed <- read.delim(bed_file, comment.char = '#', header = FALSE, fill = TRUE,
                    col.names = c('chr', 'start0', 'end', 'name'))
  bed <- bed[grepl('^chr', bed$chr), ]
  bed$STRs_ID <- sub('^[^_]+_', '', bed$name)  # name = <gene>_<STRs_ID>
  bed <- bed[bed$STRs_ID %in% variants$STRs_ID, c('STRs_ID', 'end')]
  variants$end <- bed$end[match(variants$STRs_ID, bed$STRs_ID)]
  if (anyNA(variants$end)) {
    warning('Alguns loci ausentes no BED; estimando end por motivo*copy')
    na <- is.na(variants$end)
    variants$end[na] <- variants$start0[na] + nchar(variants$motif[na]) * variants$copy[na]
  }
} else {
  cat('BED nao encontrado; estimando end por motivo*copy\n')
  variants$end <- variants$start0 + nchar(variants$motif) * variants$copy
}

# ---- leitura dos residuos normalizados (fonte dos controles) ----
norm <- readr::read_tsv(norm_file, show_col_types = FALSE)
norm$STRs_ID   <- trimws(norm$STRs_ID)
norm$sample_id <- trimws(norm$sample_id)
norm$resid     <- suppressWarnings(as.numeric(norm$allele2_residuals))

# ---- helper: acha BAM+`.bai` indexado a partir de um sample id ----
# testa `sid`, `sid.bam` e a base sem extensao; retorna o caminho do BAM ou NULL
find_bam <- function(bam_dir, sid) {
  cands <- unique(c(sid, paste0(sid, '.bam'), sub('\\.bam$', '', sid)))
  for (b in cands) {
    bfull   <- file.path(bam_dir, b)
    baifull <- sub('\\.bam$', '.bai', bfull)
    if (file.exists(bfull) && file.exists(baifull)) return(bfull)
  }
  NULL
}
# nome "limpo" (sem .bam) para o identificador do BED
clean_name <- function(sid) sub('\\.bam$', '', sid)

verify_tab <- data.frame(
  gene = character(), STRs_ID = character(),
  variant_sample = character(), variant_bam = character(),
  variant_indexed = logical(),
  control_sample = character(), control_bam = character(),
  control_indexed = logical(),
  stringsAsFactors = FALSE)

# mapeamento locus -> BAMs (caminhos absolutos) p/ carregar reads no IGV
mapping <- data.frame(
  gene = character(), STRs_ID = character(), chr = character(),
  start0 = integer(), end = integer(),
  variant_sample = character(), variant_bam = character(),
  control_sample = character(), control_bam = character(),
  stringsAsFactors = FALSE)

with_rows  <- list()
without_rows <- list()

for (i in seq_len(nrow(variants))) {
  v_strs <- variants$STRs_ID[i]
  v_chr  <- variants$chr[i]
  v_s0   <- variants$start0[i]
  v_end  <- variants$end[i]
  v_gene <- variants$gene[i]

  # ---- amostra COM a variante ----
  v_sid_raw <- as.character(variants$outlier_samples[i])
  v_bam <- find_bam(bam_dir, v_sid_raw)
  v_indexed <- !is.null(v_bam)
  v_sample  <- clean_name(v_sid_raw)

  # ---- amostra SEM a variante (controle) ----
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
    warning('Sem BAM indexado p/ amostra-com-variante: ', v_sample, ' (', v_strs, ')')
  if (!c_indexed)
    warning('Sem controle valido (BAM indexado) p/ locus: ', v_gene, ' (', v_strs, ')')

  # score: |residual| escalado (variante usa abs_res do CSV; controle usa resid)
  v_score <- if (!is.null(variants$abs_res)) {
    s <- suppressWarnings(as.numeric(variants$abs_res[i])); s <- if (is.na(s)) 0 else round(abs(s) * 100)
  } else 0
  c_score <- if (!is.na(c_sample)) {
    s <- norm[norm$STRs_ID == v_strs & clean_name(norm$sample_id) == c_sample, ]$resid
    s <- if (length(s) && !is.na(s[1])) round(abs(s[1]) * 100) else 0
  } else 0

  s0 <- max(0, v_s0 - pad)
  e0 <- v_end + pad

  if (v_indexed) {
    with_rows[[length(with_rows) + 1]] <- data.frame(
      chr = v_chr, start = s0, end = e0,
      name = paste0(v_strs, '_', v_sample),
      score = v_score, strand = '.', stringsAsFactors = FALSE)
  }
  if (c_indexed) {
    without_rows[[length(without_rows) + 1]] <- data.frame(
      chr = v_chr, start = s0, end = e0,
      name = paste0(v_strs, '_', c_sample),
      score = c_score, strand = '.', stringsAsFactors = FALSE)
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
    cat('Nenhuma linha valida p/', out, '- arquivo NAO escrito.\n')
    return(0)
  }
  df <- do.call(rbind, rows)
  df <- df[order(df$chr, df$start), ]
  write.table(df, out, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat('Escrevi', nrow(df), 'loci ->', out, '\n')
  nrow(df)
}

n_with    <- write_bed(with_rows, out_with)
n_without <- write_bed(without_rows, out_without)

out_map <- 'str_samples_bams.tsv'
write.table(mapping, out_map, sep = '\t', row.names = FALSE, quote = FALSE)
cat('Escrevi', nrow(mapping), 'loci ->', out_map, '\n')

cat('\n=== Resumo (verify_tab) ===\n')
print(verify_tab)

cat('\nNo IGV: File > Load from File ->', out_with, 'e', out_without, '\n')
cat('Dica: sobreponha aos bigWigs/bigBeds de external_tracks/.\n')
