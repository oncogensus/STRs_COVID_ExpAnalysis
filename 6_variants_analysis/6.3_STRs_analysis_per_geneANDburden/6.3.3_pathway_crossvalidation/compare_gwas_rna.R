# compare_gwas_rna.R
# ---------------------------------------------------------------------------
# Comparacao entre as estrategias GWAS-filtrado e RNA-seq:
#   1) outliers  : STRs/genes com STR-outlier por estrategia + uniao
#   2) sobreposicao por STR x paciente (grupo caso vs controle) - DESCRITIVO
#      (sem testes), incluindo overlap do tamanho do maior alelo entre grupos
#
# Entradas (padroes cluster):
#   .../6.3.1.2.../results/covid_suggestive_genes_with_outlier_STRs.tsv (P1 GWAS)
#   .../6.3.2.1_RNA_matrix/results/rna_outlier_genes.tsv                (outliers RNA)
#   <repo>/samples/STRs_analysis_dataset.tsv                            (STR x paciente)
#
# Saidas (--out-dir, padrao results_gwas_rna_comparison/):
#   strategy_outlier_sets.tsv / outlier_genes_union.tsv
#   patient_str.tsv / per_str_case_control.tsv
#
# Uso:
#   Rscript compare_gwas_rna.R [--repo <dir>] [--out-dir <dir>]
#                              [--p1-file ...] [--rna-outliers ...]
#                              [--catalog ...]
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
  file.path(TOP, "6.3.1_GWAS_analysis/6.3.1.2_dbscan_subset_GWAS/results/covid_suggestive_genes_with_outlier_STRs.tsv"))
rna_out       <- get_opt(cmd_args, "--rna-outliers",
  file.path(TOP, "6.3.2_RNA_data_analysis/6.3.2.1_RNA_matrix/results/rna_outlier_genes.tsv"))
catalog_file  <- get_opt(cmd_args, "--catalog",
  file.path(REPO, "samples/STRs_analysis_dataset.tsv"))
out_dir       <- get_opt(cmd_args, "--out-dir",
  file.path(TOP, "6.3.3_pathway_crossvalidation/results_gwas_rna_comparison"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (p in c(p1_file, rna_out, catalog_file))
  if (!file.exists(p)) stop("arquivo ausente: ", p)

num <- function(x) suppressWarnings(as.numeric(x))

cat("=== Comparacao GWAS-filtrado x RNA-seq ===\n")
cat("p1_file      :", p1_file, "\n")
cat("rna_outliers :", rna_out, "\n")
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
## 2. CONJUNTOS DE STR POR ESTRATEGIA
## ---------------------------------------------------------------------------
union_strs <- unique(c(gwas_strs, rna_strs))

marker <- function(ids) data.table(strs_id = unique(ids))

sets <- unique(rbindlist(list(
  gwas_meta[, .(strs_id, gene)],
  rna_meta[, .(strs_id, gene)]
)))
sets[, in_gwas_sig := 0L][marker(gwas_strs), in_gwas_sig := 1L, on = "strs_id"]
sets[, in_rna := 0L][marker(rna_strs), in_rna := 1L, on = "strs_id"]

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
cat(sprintf("Uniao STRs: %d | in_gwas_sig=%d, in_rna=%d\n",
            uniqueN(sets$strs_id),
            length(gwas_strs), length(rna_strs)))

## ---------------------------------------------------------------------------
## 3. UNIAO DE GENES POR ESTRATEGIA
## ---------------------------------------------------------------------------
pairs <- unique(rbindlist(list(
  sets[in_gwas_sig == 1L, .(gene, strs_id, src = "gwas_sig")],
  sets[in_rna == 1L, .(gene, strs_id, src = "rna")]
), use.names = TRUE))

genes_union <- dcast(pairs, gene ~ src, value.var = "strs_id",
                     fun.aggregate = length)
src_new <- c(gwas_sig = "n_strs_gwas_sig", rna = "n_strs_rna")
rn <- intersect(names(src_new), names(genes_union))
setnames(genes_union, rn, src_new[rn])
for (cc in c("n_strs_gwas_sig", "n_strs_rna")) {
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
setorder(genes_union, -in_gwas_sig, -in_rna, gene)

fwrite(genes_union, file.path(out_dir, "outlier_genes_union.tsv"), sep = "\t")
cat(sprintf("Genes na uniao: %d | com outlier GWAS sig=%d, RNA=%d\n",
            nrow(genes_union),
            sum(genes_union$in_gwas_sig == 1L),
            sum(genes_union$in_rna == 1L)))

## ---------------------------------------------------------------------------
## 4. TABELA LONGA STR x PACIENTE (do catalogo da coorte)
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

setorder(pat, strs_id, sample_id)
fwrite(pat, file.path(out_dir, "patient_str.tsv"), sep = "\t")
cat(sprintf("Tabela longa STR x paciente: %d linhas | %d pacientes | %d grupos\n",
            nrow(pat), uniqueN(pat$sample_id),
            uniqueN(pat$group[!is.na(pat$group)])))

## ---------------------------------------------------------------------------
## 5. RESUMO DESCRITIVO POR STR (caso x controle, SEM testes)
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
## 6. OVERLAP DO MAIOR ALELO POR GENE (agregado dos STRs)
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
setorder(genes_union, -in_gwas_sig, -in_rna, gene)
fwrite(genes_union, file.path(out_dir, "outlier_genes_union.tsv"), sep = "\t")

## ---------------------------------------------------------------------------
## 7. RESUMO FINAL
## ---------------------------------------------------------------------------
cat("\n=== RESULTADO DESCRITIVO POR STR (caso x controle) ===\n")
cat(sprintf("  STRs na uniao: %d\n", uniqueN(per_str$strs_id)))
cat(sprintf("  overlap maior alelo 'sim': %d | 'nao': %d | 'sem_dados': %d\n",
            sum(per_str$overlap_maior_alealo_grupos == "sim"),
            sum(per_str$overlap_maior_alealo_grupos == "nao"),
            sum(per_str$overlap_maior_alealo_grupos == "sem_dados")))
for (src in c("in_gwas_sig", "in_rna")) {
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
