# compare_gwas_rna.R
# ---------------------------------------------------------------------------
# Script unico de comparacao entre as estrategias GWAS-filtrado e RNA-seq:
#   1) outliers  : STRs/genes com STR-outlier por estrategia + uniao
#   2) burden    : hits SKAT por gene (uncorrected p<0.05, corrected BH q<0.05)
#   3) sobreposicao por STR x paciente (grupo caso vs controle) - DESCRITIVO
#      (sem testes), incluindo overlap do tamanho do maior alelo entre grupos
#
# Entradas (padroes cluster):
#   .../6.4.1.2.../results/covid_suggestive_genes_with_outlier_STRs.tsv (P1 GWAS)
#   .../6.4.2.2_RNA_matrix/results/rna_outlier_genes.tsv                (outliers RNA)
#   .../6.4.3_burden_test/results_gwas_burden/skat_per_gene.tsv         (burden GWAS)
#   .../6.4.3_burden_test/results_rna/skat_per_gene.tsv                 (burden RNA)
#   <repo>/samples/STRs_analysis_dataset.tsv                            (STR x paciente)
#
# Saidas (--out-dir, padrao results_gwas_rna_comparison/):
#   strategy_outlier_sets.tsv / outlier_genes_union.tsv
#   burden_hits_union_{uncorrected,corrected}.tsv / burden_hits_overlap_{...}.tsv
#   outlier_x_burden_by_gene.tsv
#   patient_str.tsv / per_str_case_control.tsv
#
# Uso:
#   Rscript compare_gwas_rna.R [--repo <dir>] [--out-dir <dir>]
#                              [--p1-file ...] [--rna-outliers ...]
#                              [--gwas-skat ...] [--rna-skat ...] [--catalog ...]
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table) })

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

cmd_args <- commandArgs(trailingOnly = TRUE)
REPO <- get_opt(cmd_args, "--repo",
                "/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis")
TOP <- file.path(REPO, "6_variants_analysis", "6.3_STRs_analysis_per_geneANDburden")

p1_file       <- get_opt(cmd_args, "--p1-file",
  file.path(TOP, "6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv"))
rna_out       <- get_opt(cmd_args, "--rna-outliers",
  file.path(TOP, "6.4.2_dbscan_subset_RNA/6.4.2.2_RNA_matrix/results/rna_outlier_genes.tsv"))
gwas_skat_file <- get_opt(cmd_args, "--gwas-skat",
  file.path(TOP, "6.4.3_burden_test/results_gwas_burden/skat_per_gene.tsv"))
rna_skat_file  <- get_opt(cmd_args, "--rna-skat",
  file.path(TOP, "6.4.3_burden_test/results_rna/skat_per_gene.tsv"))
catalog_file  <- get_opt(cmd_args, "--catalog",
  file.path(REPO, "samples/STRs_analysis_dataset.tsv"))
out_dir       <- get_opt(cmd_args, "--out-dir",
  file.path(TOP, "6.4.4_pathway_crossvalidation/results_gwas_rna_comparison"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (p in c(p1_file, rna_out, gwas_skat_file, rna_skat_file, catalog_file))
  if (!file.exists(p)) stop("arquivo ausente: ", p)

num <- function(x) suppressWarnings(as.numeric(x))

cat("=== Comparacao GWAS-filtrado x RNA-seq ===\n")
cat("p1_file      :", p1_file, "\n")
cat("rna_outliers :", rna_out, "\n")
cat("gwas_skat    :", gwas_skat_file, "\n")
cat("rna_skat     :", rna_skat_file, "\n")
cat("catalog      :", catalog_file, "\n")
cat("out_dir      :", out_dir, "\n")

## ---------------------------------------------------------------------------
## 1. OUTLIERS POR ESTRATEGIA
## ---------------------------------------------------------------------------
p1 <- fread(p1_file, header = TRUE, sep = "\t")
p1[, sn_out := num(subset_n_outliers)]
gwas_out <- p1[gwas_significance == "significant" & is.finite(sn_out) & sn_out >= 1]

rna_outl <- fread(rna_out, header = TRUE, sep = "\t")

gwas_strs <- unique(gwas_out$strs_id)
rna_strs  <- unique(rna_outl$strs_id)
cat(sprintf("Outliers GWAS (significativos): %d STRs | %d genes\n",
            length(gwas_strs), uniqueN(gwas_out$gene)))
cat(sprintf("Outliers RNA:                    %d STRs | %d genes\n",
            length(rna_strs), uniqueN(rna_outl$gene)))

gwas_meta <- unique(gwas_out[, .(strs_id, gene, region, repeat_unit,
                                 gwas_p = num(gwas_p), n_out_gwas = sn_out)])
rna_meta  <- unique(rna_outl[, .(strs_id, gene, region, repeat_unit,
                                 n_out_rna = num(n_outliers_dbscan_global))])

gwas_samp <- unique(gwas_out[, .(strs_id, os = subset_outlier_samples)])
rna_samp  <- unique(rna_outl[, .(strs_id, os = outlier_samples_dbscan_global)])

## ---------------------------------------------------------------------------
## 2. BURDEN: STRs e genes testados / hits
## ---------------------------------------------------------------------------
burden_read <- function(path) {
  d <- fread(path, header = TRUE, sep = "\t")
  if (!"gene" %in% names(d))  stop("coluna 'gene' ausente em: ", path)
  if (!"str_ids" %in% names(d)) stop("coluna 'str_ids' ausente em: ", path)
  d
}
gwas_skat <- burden_read(gwas_skat_file)
rna_skat  <- burden_read(rna_skat_file)

strs_of_genes <- function(d) {
  ids <- trimws(unlist(strsplit(as.character(d$str_ids), ";")))
  unique(ids[nzchar(ids)])
}
burden_gwas_strs <- strs_of_genes(gwas_skat)
burden_rna_strs  <- strs_of_genes(rna_skat)
cat(sprintf("Burden GWAS: %d genes testados | %d STRs | %d hits unc / %d cor\n",
            uniqueN(gwas_skat$gene), length(burden_gwas_strs),
            sum(!is.na(gwas_skat$p_value) & gwas_skat$p_value < 0.05, na.rm = TRUE),
            sum(!is.na(gwas_skat$q_value) & gwas_skat$q_value < 0.05, na.rm = TRUE)))
cat(sprintf("Burden RNA:  %d genes testados | %d STRs | %d hits unc / %d cor\n",
            uniqueN(rna_skat$gene), length(burden_rna_strs),
            sum(!is.na(rna_skat$p_value) & rna_skat$p_value < 0.05, na.rm = TRUE),
            sum(!is.na(rna_skat$q_value) & rna_skat$q_value < 0.05, na.rm = TRUE)))

skat_gene_map <- function(d) {
  out <- lapply(seq_len(nrow(d)), function(i) {
    ids <- trimws(unlist(strsplit(as.character(d$str_ids[i]), ";")))
    ids <- ids[nzchar(ids)]
    gn <- as.character(d$gene_name[i])
    if (!length(ids) || is.na(gn)) return(NULL)
    data.table(strs_id = ids, gene = gn)
  })
  rbindlist(Filter(Negate(is.null), out))
}
skat_map <- unique(rbindlist(list(skat_gene_map(gwas_skat),
                                  skat_gene_map(rna_skat))))

## ---------------------------------------------------------------------------
## 3. CONJUNTOS DE STR POR ESTRATEGIA
## ---------------------------------------------------------------------------
union_strs <- unique(c(gwas_strs, rna_strs, burden_gwas_strs, burden_rna_strs))

marker <- function(ids) data.table(strs_id = unique(ids))

sets <- unique(rbindlist(list(
  gwas_meta[, .(strs_id, gene)],
  rna_meta[, .(strs_id, gene)],
  skat_map[, .(strs_id, gene)]
)))
sets[, in_gwas_sig := 0L][marker(gwas_strs), in_gwas_sig := 1L, on = "strs_id"]
sets[, in_rna := 0L][marker(rna_strs), in_rna := 1L, on = "strs_id"]
sets[, in_burden_gwas := 0L][marker(burden_gwas_strs), in_burden_gwas := 1L, on = "strs_id"]
sets[, in_burden_rna := 0L][marker(burden_rna_strs), in_burden_rna := 1L, on = "strs_id"]

str_info <- unique(rbindlist(list(
  gwas_meta[, .(strs_id, region, repeat_unit, gwas_p, n_out_gwas)],
  rna_meta[, .(strs_id, region, repeat_unit, n_out_rna)]
), use.names = TRUE, fill = TRUE))
str_info[, gwas_p     := if (all(is.na(gwas_p))) NA_real_ else min(gwas_p[!is.na(gwas_p)]), by = strs_id]
str_info[, n_out_gwas := if (all(is.na(n_out_gwas))) NA_real_ else max(n_out_gwas[!is.na(n_out_gwas)]), by = strs_id]
str_info[, n_out_rna  := if (all(is.na(n_out_rna)))  NA_real_ else max(n_out_rna[!is.na(n_out_rna)]),  by = strs_id]
str_info[, region      := region[1],      by = strs_id]
str_info[, repeat_unit := repeat_unit[1], by = strs_id]
str_info <- unique(str_info, by = "strs_id")

sets <- merge(sets, str_info, by = "strs_id", all.x = TRUE)
setorder(sets, gene, strs_id)

fwrite(sets, file.path(out_dir, "strategy_outlier_sets.tsv"), sep = "\t")
cat(sprintf("Uniao STRs: %d | in_gwas_sig=%d, in_rna=%d, in_burden_gwas=%d, in_burden_rna=%d\n",
            uniqueN(sets$strs_id),
            length(gwas_strs), length(rna_strs),
            length(burden_gwas_strs), length(burden_rna_strs)))

## ---------------------------------------------------------------------------
## 4. UNIAO DE GENES POR ESTRATEGIA
## ---------------------------------------------------------------------------
pairs <- unique(rbindlist(list(
  sets[in_gwas_sig == 1L, .(gene, strs_id, src = "gwas_sig")],
  sets[in_rna == 1L, .(gene, strs_id, src = "rna")],
  sets[in_burden_gwas == 1L, .(gene, strs_id, src = "burden_gwas")],
  sets[in_burden_rna == 1L, .(gene, strs_id, src = "burden_rna")]
), use.names = TRUE))

genes_union <- dcast(pairs, gene ~ src, value.var = "strs_id",
                     fun.aggregate = length)
src_new <- c(gwas_sig = "n_strs_gwas_sig", rna = "n_strs_rna",
             burden_gwas = "n_strs_burden_gwas", burden_rna = "n_strs_burden_rna")
rn <- intersect(names(src_new), names(genes_union))
setnames(genes_union, rn, src_new[rn])
for (cc in c("n_strs_gwas_sig", "n_strs_rna", "n_strs_burden_gwas", "n_strs_burden_rna")) {
  if (is.null(genes_union[[cc]])) genes_union[[cc]] <- 0L
  genes_union[[cc]][is.na(genes_union[[cc]])] <- 0L
  genes_union[[cc]] <- as.integer(genes_union[[cc]])
}

gwas_p_gene <- gwas_out[, .(gwas_p_min = min(num(gwas_p), na.rm = TRUE)), by = gene]
genes_union <- merge(genes_union, gwas_p_gene, by = "gene", all.x = TRUE)
rna_studies <- rna_outl[, .(rna_gse = paste(sort(unique(gse)), collapse = ";")), by = gene]
genes_union <- merge(genes_union, rna_studies, by = "gene", all.x = TRUE)

genes_union[, in_gwas_sig := as.integer(n_strs_gwas_sig > 0)]
genes_union[, in_rna := as.integer(n_strs_rna > 0)]
genes_union[, in_burden_gwas := as.integer(n_strs_burden_gwas > 0)]
genes_union[, in_burden_rna := as.integer(n_strs_burden_rna > 0)]
setorder(genes_union, -in_gwas_sig, -in_rna, -in_burden_gwas, gene)

fwrite(genes_union, file.path(out_dir, "outlier_genes_union.tsv"), sep = "\t")
cat(sprintf("Genes na uniao: %d | com outlier GWAS sig=%d, RNA=%d\n",
            nrow(genes_union),
            sum(genes_union$in_gwas_sig == 1L),
            sum(genes_union$in_rna == 1L)))

## ---------------------------------------------------------------------------
## 5. SOBREPOSICAO DE HITS DO BURDEN (union / overlap, uncorrected e corrected)
## ---------------------------------------------------------------------------
build_union <- function(g, r, hit_g_col, hit_r_col) {
  keep_g <- c("gene", "gene_name", "n_variants", "str_ids",
              "p_value", "q_value", "note", hit_g_col)
  keep_r <- c("gene", "gene_name", "n_variants", "str_ids",
              "p_value", "q_value", "note", hit_r_col)
  g <- g[, keep_g, drop = FALSE]
  r <- r[, keep_r, drop = FALSE]

  names(g)[names(g) == "gene_name"]  <- "gwas_gene_name"
  names(g)[names(g) == "n_variants"] <- "gwas_n_variants"
  names(g)[names(g) == "str_ids"]    <- "gwas_str_ids"
  names(g)[names(g) == "p_value"]    <- "gwas_p_value"
  names(g)[names(g) == "q_value"]    <- "gwas_q_value"
  names(g)[names(g) == "note"]       <- "gwas_note"
  names(r)[names(r) == "gene_name"]  <- "rna_gene_name"
  names(r)[names(r) == "n_variants"] <- "rna_n_variants"
  names(r)[names(r) == "str_ids"]    <- "rna_str_ids"
  names(r)[names(r) == "p_value"]    <- "rna_p_value"
  names(r)[names(r) == "q_value"]    <- "rna_q_value"
  names(r)[names(r) == "note"]       <- "rna_note"

  both <- merge(g, r, by = "gene", all = TRUE)
  both$gene_name <- ifelse(is.na(both$gwas_gene_name),
                           both$rna_gene_name, both$gwas_gene_name)
  both$gwas_gene_name <- NULL
  both$rna_gene_name <- NULL
  both[[hit_g_col]][is.na(both[[hit_g_col]])] <- 0L
  both[[hit_r_col]][is.na(both[[hit_r_col]])] <- 0L

  sort_p_r <- ifelse(is.na(both$rna_p_value), 1, both$rna_p_value)
  sort_p_g <- ifelse(is.na(both$gwas_p_value), 1, both$gwas_p_value)
  both$priority <- both[[hit_g_col]] + both[[hit_r_col]]
  both <- both[order(-both$priority, sort_p_r, sort_p_g), ]
  both$priority <- NULL

  both[, c("gene", "gene_name", hit_g_col, hit_r_col,
           "gwas_n_variants", "rna_n_variants",
           "gwas_p_value", "rna_p_value", "gwas_q_value", "rna_q_value",
           "gwas_str_ids", "rna_str_ids"), drop = FALSE]
}

mk_hits <- function(d) {
  list(unc = d[!is.na(p_value) & p_value < 0.05],
       cor = d[!is.na(q_value) & q_value < 0.05])
}
hits_df <- function(dt) as.data.frame(dt)

g_unc <- hits_df(mk_hits(gwas_skat)$unc)
g_cor <- hits_df(mk_hits(gwas_skat)$cor)
r_unc <- hits_df(mk_hits(rna_skat)$unc)
r_cor <- hits_df(mk_hits(rna_skat)$cor)
g_unc$hit_gwas <- rep(1L, nrow(g_unc))
g_cor$hit_gwas <- rep(1L, nrow(g_cor))
r_unc$hit_rna  <- rep(1L, nrow(r_unc))
r_cor$hit_rna  <- rep(1L, nrow(r_cor))

burden_unc <- build_union(g_unc, r_unc, "hit_gwas", "hit_rna")
burden_cor <- build_union(g_cor, r_cor, "hit_gwas", "hit_rna")

ov_unc <- burden_unc[burden_unc$hit_gwas == 1L & burden_unc$hit_rna == 1L, , drop = FALSE]
ov_cor <- burden_cor[burden_cor$hit_gwas == 1L & burden_cor$hit_rna == 1L, , drop = FALSE]

fwrite(burden_unc, file.path(out_dir, "burden_hits_union_uncorrected.tsv"), sep = "\t")
fwrite(burden_cor, file.path(out_dir, "burden_hits_union_corrected.tsv"), sep = "\t")
fwrite(ov_unc, file.path(out_dir, "burden_hits_overlap_uncorrected.tsv"), sep = "\t")
fwrite(ov_cor, file.path(out_dir, "burden_hits_overlap_corrected.tsv"), sep = "\t")

cat("Sobreposicao burden (uncorrected p<0.05): GWAS=", sum(burden_unc$hit_gwas == 1L),
    " RNA=", sum(burden_unc$hit_rna == 1L), " overlap=", nrow(ov_unc), "\n", sep = "")
cat("Sobreposicao burden (corrected q<0.05):   GWAS=", sum(burden_cor$hit_gwas == 1L),
    " RNA=", sum(burden_cor$hit_rna == 1L), " overlap=", nrow(ov_cor), "\n", sep = "")

## ---------------------------------------------------------------------------
## 6. MATRIZ 2x2 OUTLIER x BURDEN POR GENE
## ---------------------------------------------------------------------------
bur_hits_gene <- function(d) {
  list(unc = unique(as.character(d[!is.na(p_value) & p_value < 0.05]$gene_name)),
       cor = unique(as.character(d[!is.na(q_value) & q_value < 0.05]$gene_name)))
}
bh_g <- bur_hits_gene(gwas_skat)
bh_r <- bur_hits_gene(rna_skat)

out_g <- unique(gwas_out$gene)
out_r <- unique(rna_outl$gene)

oxb_genes <- unique(c(out_g, out_r, bh_g$unc, bh_r$unc, bh_g$cor, bh_r$cor))
oxb <- data.table(gene   = oxb_genes,
                  out_gwas_sig = as.integer(oxb_genes %in% out_g),
                  out_rna       = as.integer(oxb_genes %in% out_r),
                  bur_gwas_unc  = as.integer(oxb_genes %in% bh_g$unc),
                  bur_rna_unc   = as.integer(oxb_genes %in% bh_r$unc),
                  bur_gwas_cor  = as.integer(oxb_genes %in% bh_g$cor),
                  bur_rna_cor   = as.integer(oxb_genes %in% bh_r$cor))
setorder(oxb, -bur_gwas_unc, -bur_rna_unc, -out_gwas_sig, -out_rna, gene)
fwrite(oxb, file.path(out_dir, "outlier_x_burden_by_gene.tsv"), sep = "\t")

cat("\n=== 2x2 Outlier x Burden (por gene) ===\n")
cat(sprintf("  out_gwas_sig & bur_gwas_unc: %d | out_gwas_sig & bur_rna_unc: %d\n",
            sum(oxb$out_gwas_sig & oxb$bur_gwas_unc),
            sum(oxb$out_gwas_sig & oxb$bur_rna_unc)))
cat(sprintf("  out_rna      & bur_rna_unc : %d | out_rna      & bur_gwas_unc: %d\n",
            sum(oxb$out_rna & oxb$bur_rna_unc),
            sum(oxb$out_rna & oxb$bur_gwas_unc)))
cat(sprintf("  com overlap corrigido (q<0.05): out_gwas_sig=%d out_rna=%d\n",
            sum(oxb$out_gwas_sig & oxb$bur_gwas_cor),
            sum(oxb$out_rna & oxb$bur_rna_cor)))

## ---------------------------------------------------------------------------
## 7. TABELA LONGA STR x PACIENTE (do catalogo da coorte)
## ---------------------------------------------------------------------------
cat_cols <- c("STRs_ID", "sample_id", "group", "allele1_est", "allele2_est",
              "chrom", "start", "end", "region", "repeat_unit")
catc <- fread(catalog_file, header = TRUE, sep = "\t", select = cat_cols)
cat(sprintf("Catalogo carregado: %d linhas\n", nrow(catc)))

union_sid <- marker(union_strs)
pat <- catc[union_sid, on = c("STRs_ID" = "strs_id"), nomatch = 0L]
setnames(pat, "STRs_ID", "strs_id")

pat[, sample_id := as.character(sample_id)]
pat[, group := tolower(trimws(as.character(group)))]
pat[, allele1_est := num(allele1_est)]
pat[, allele2_est := num(allele2_est)]
pat[, maior_alelo := ifelse(is.na(allele1_est) & is.na(allele2_est), NA_real_,
                            pmax(allele1_est, allele2_est, na.rm = TRUE))]
pat <- unique(pat, by = c("strs_id", "sample_id"))
pat[, chrom := NULL][, start := NULL][, end := NULL]

pat[, in_gwas_sig := 0L][marker(gwas_strs), in_gwas_sig := 1L, on = "strs_id"]
pat[, in_rna := 0L][marker(rna_strs), in_rna := 1L, on = "strs_id"]
pat[, in_burden_gwas := 0L][marker(burden_gwas_strs), in_burden_gwas := 1L, on = "strs_id"]
pat[, in_burden_rna := 0L][marker(burden_rna_strs), in_burden_rna := 1L, on = "strs_id"]

setorder(pat, strs_id, sample_id)
fwrite(pat, file.path(out_dir, "patient_str.tsv"), sep = "\t")
cat(sprintf("Tabela longa STR x paciente: %d linhas | %d pacientes | %d grupos\n",
            nrow(pat), uniqueN(pat$sample_id),
            uniqueN(pat$group[!is.na(pat$group)])))

## ---------------------------------------------------------------------------
## 8. RESUMO DESCRITIVO POR STR (caso x controle, SEM testes)
## ---------------------------------------------------------------------------
parse_ids <- function(x) {
  if (is.na(x) || is.null(x)) return(character(0))
  ids <- trimws(unlist(strsplit(as.character(x), "[;,]")[[1]]))
  ids[nzchar(ids)]
}

sum_stats <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) {
    return(list(n = 0L, mean = NA_real_, median = NA_real_,
                min = NA_real_, max = NA_real_, sd = NA_real_))
  }
  list(n = length(v), mean = mean(v), median = median(v),
       min = min(v), max = max(v), sd = if (length(v) > 1) sd(v) else NA_real_)
}

count_out_grp <- function(ids, pat_s) {
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (!length(ids)) return(c(case = 0L, control = 0L))
  h <- pat_s[sample_id %in% ids]
  c(case = sum(h$group == "case"), control = sum(h$group == "control"))
}

gwas_os_key <- if (nrow(gwas_samp)) unique(gwas_samp, by = "strs_id") else data.table()
rna_os_key  <- if (nrow(rna_samp))  unique(rna_samp, by = "strs_id")  else data.table()

per_str_list <- lapply(unique(pat$strs_id), function(sid) {
  info <- sets[strs_id == sid, ][1]
  if (!nrow(info)) return(NULL)
  p1l <- pat[strs_id == sid]
  sc <- sum_stats(p1l[group == "case"]$maior_alelo)
  sn <- sum_stats(p1l[group == "control"]$maior_alelo)

  os_g <- if (nrow(gwas_os_key)) gwas_os_key[strs_id == sid]$os[1] else NA_character_
  os_r <- if (nrow(rna_os_key))  rna_os_key[strs_id == sid]$os[1]  else NA_character_
  og <- count_out_grp(parse_ids(os_g), p1l)
  orr <- count_out_grp(parse_ids(os_r), p1l)

  ov <- "sem_dados"
  if (sc$n && sn$n) ov <- if (sc$min <= sn$max && sn$min <= sc$max) "sim" else "nao"

  data.table(
    strs_id = sid,
    gene = info$gene,
    region = info$region,
    repeat_unit = info$repeat_unit,
    in_gwas_sig = info$in_gwas_sig, in_rna = info$in_rna,
    in_burden_gwas = info$in_burden_gwas, in_burden_rna = info$in_burden_rna,
    n_case = sc$n, n_control = sn$n,
    n_out_gwas_case = og["case"], n_out_gwas_control = og["control"],
    n_out_rna_case = orr["case"], n_out_rna_control = orr["control"],
    mean_case = sc$mean, median_case = sc$median, min_case = sc$min,
    max_case = sc$max, sd_case = sc$sd,
    mean_control = sn$mean, median_control = sn$median, min_control = sn$min,
    max_control = sn$max, sd_control = sn$sd,
    overlap_maior_alealo_grupos = ov)
})

per_str <- rbindlist(per_str_list)
setorder(per_str, strs_id)
fwrite(per_str, file.path(out_dir, "per_str_case_control.tsv"), sep = "\t")

## ---------------------------------------------------------------------------
## 9. OVERLAP DO MAIOR ALELO POR GENE (agregado dos STRs)
## ---------------------------------------------------------------------------
gene_overlap <- per_str[, .(flags = paste(overlap_maior_alealo_grupos, collapse = ";")),
                        by = gene]
gene_overlap[, overlap_maior_alealo_grupos := ifelse(
  !grepl("sim|nao", flags), "sem_dados",
  ifelse(vapply(strsplit(flags, ";"), function(z) "nao" %in% z, logical(1L)),
         "nao", "sim"))]
genes_union <- merge(genes_union,
                     gene_overlap[, .(gene, overlap_maior_alealo_grupos)],
                     by = "gene", all.x = TRUE)
setorder(genes_union, -in_gwas_sig, -in_rna, -in_burden_gwas, gene)
fwrite(genes_union, file.path(out_dir, "outlier_genes_union.tsv"), sep = "\t")

## ---------------------------------------------------------------------------
## 10. RESUMO FINAL
## ---------------------------------------------------------------------------
cat("\n=== RESULTADO DESCRITIVO POR STR (caso x controle) ===\n")
cat(sprintf("  STRs na uniao: %d\n", uniqueN(per_str$strs_id)))
cat(sprintf("  overlap maior alelo 'sim': %d | 'nao': %d | 'sem_dados': %d\n",
            sum(per_str$overlap_maior_alealo_grupos == "sim"),
            sum(per_str$overlap_maior_alealo_grupos == "nao"),
            sum(per_str$overlap_maior_alealo_grupos == "sem_dados")))
for (src in c("in_gwas_sig", "in_rna", "in_burden_gwas", "in_burden_rna")) {
  sub <- per_str[get(src) == 1L]
  if (!nrow(sub)) next
  wm <- function(m, n) if (sum(n) > 0) sum(m * n) / sum(n) else NaN
  cat(sprintf("  [%s] STRs=%d | overlap sim=%d nao=%d sem_dados=%d | caso n=%d (maior alelo med=%.1f) control n=%d (maior alelo med=%.1f)\n",
              src, nrow(sub),
              sum(sub$overlap_maior_alealo_grupos == "sim"),
              sum(sub$overlap_maior_alealo_grupos == "nao"),
              sum(sub$overlap_maior_alealo_grupos == "sem_dados"),
              sum(sub$n_case), wm(sub$mean_case, sub$n_case),
              sum(sub$n_control), wm(sub$mean_control, sub$n_control)))
}

cat("\n=== Saidas em:", out_dir, "===\n")
cat(paste(sort(list.files(out_dir)), collapse = "\n  "), "\n")
cat("\n=== FIM compare_gwas_rna.R ===\n")