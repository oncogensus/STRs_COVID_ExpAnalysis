#!/usr/bin/env Rscript
# str_regions_to_bed.R
# Extrai as regioes das STRs avaliadas como BED para visualizacao no IGV.
# Standalone: NAO depende de trackViewer / GenomicRanges.
# Rodar a partir de 7.4_scRNA-seq_analysis/ (onde STR_variants_UCSC_track.bed esta).

variants_file <- '../7.3_STRs_filter/results/STR_vs_scRNA_overlap_unified.csv'
bed_file      <- 'STR_variants_UCSC_track.bed'
out_exact    <- 'str_regions.bed'           # regiao exata da STR (0-based start)
out_padded   <- 'str_regions.padded20kb.bed' # regiao +-20kb (contexto p/ IGV)
pad          <- 20000

stopifnot(file.exists(variants_file))

# CSV (readr se disponivel, senao utils)
if (requireNamespace('readr', quietly = TRUE)) {
  ov <- readr::read_csv(variants_file, show_col_types = FALSE)
} else {
  ov <- read.csv(variants_file, stringsAsFactors = FALSE)
}

# uma linha por locus de STR (remove linhas duplicadas de cell_type)
variants <- ov[!duplicated(ov$STRs_ID), ]

# STRs_ID codifica chr:pos:motif:copy  (pos = 1-based start)
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

# score = |residual absoluto| escalado (maior = mais extremo); 0 se ausente
score <- if (!is.null(variants$abs_res)) round(abs(as.numeric(variants$abs_res)) * 100) else 0
score[is.na(score)] <- 0
# nome: gene|STRs_ID (cai no STRs_ID se gene vazio)
name <- ifelse(is.na(variants$gene) | variants$gene == '',
               variants$STRs_ID, paste0(variants$gene, '|', variants$STRs_ID))

bed6 <- data.frame(chr = variants$chr,
                   start = variants$start0,        # 0-based (BED)
                   end = variants$end,             # 1-based-style (BED)
                   name = name, score = score, strand = '.',
                   stringsAsFactors = FALSE)
bed6 <- bed6[order(bed6$chr, bed6$start), ]
write.table(bed6, out_exact, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
cat('Escrevi', nrow(bed6), 'loci de STR ->', out_exact, '\n')

# versao com +-pad kb de contexto para navegar no IGV
pad_df <- data.frame(chr = variants$chr,
                     start = pmax(0, variants$start0 - pad),
                     end = variants$end + pad,
                     name = name, score = score, strand = '.',
                     stringsAsFactors = FALSE)
pad_df <- pad_df[order(pad_df$chr, pad_df$start), ]
write.table(pad_df, out_padded, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
cat('Escrevi regioes com contexto (+-', pad / 1e3, 'kb) ->', out_padded, '\n')

cat('\nExemplo (cabecalho):\n'); print(head(bed6))
cat('\nNo IGV: File > Load from File ->', out_exact, '(ou', out_padded, ')\n')
cat('Dica: os bigWigs/bigBeds de external_tracks tambem abrem direto no IGV (File > Load from File).\n')
