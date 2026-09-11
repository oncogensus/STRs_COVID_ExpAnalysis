# compare_gwas_rna.R
# ---------------------------------------------------------------------------
# Analise descritiva dos outliers RNA-seq: STRs/genes com outlier DBSCAN
# global, sobreposicao por STR x paciente (grupo caso vs controle) - DESCRITIVO
# (sem testes), incluindo overlap do tamanho do maior alelo entre grupos.
#
# Entradas (padroes cluster):
#   .../6.3.2.1_RNA_matrix/results/rna_outlier_genes.tsv  (outliers RNA)
#   <repo>/samples/STRs_analysis_dataset.tsv              (STR x paciente)
#
# Saidas (--out-dir, padrao results_gwas_rna_comparison/):
#   rna_outlier_sets.tsv / rna_genes_summary.tsv
#   patient_str.tsv / per_str_case_control.tsv
#
# Uso:
#   Rscript compare_gwas_rna.R [--repo <dir>] [--out-dir <dir>]
#                              [--rna-outliers ...] [--catalog ...]
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

rna_out       <- get_opt(cmd_args, "--rna-outliers",
  file.path(TOP, "6.3.2_RNA_data_analysis/6.3.2.1_RNA_matrix/results/rna_outlier_genes.tsv"))
catalog_file  <- get_opt(cmd_args, "--catalog",
  file.path(REPO, "samples/STRs_analysis_dataset.tsv"))
out_dir       <- get_opt(cmd_args, "--out-dir",
  file.path(TOP, "6.3.3_pathway_crossvalidation/results_gwas_rna_comparison"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (p in c(rna_out, catalog_file))
  if (!file.exists(p)) stop("arquivo ausente: ", p)

num <- function(x) suppressWarnings(as.numeric(x))

cat("=== Analise descritiva outliers RNA-seq ===\n")
cat("rna_outliers :", rna_out, "\n")
cat("catalog      :", catalog_file, "\n")
cat("out_dir      :", out_dir, "\n")

## ---------------------------------------------------------------------------
## 1. OUTLIERS RNA-SEQ
## ---------------------------------------------------------------------------
rna_outl <- fread(rna_out, header = TRUE, sep = "\t")

rna_strs  <- unique(rna_outl$strs_id)
cat(sprintf("Outliers RNA: %d STRs | %d genes\n",
            length(rna_strs), uniqueN(rna_outl$gene)))

rna_meta <- unique(rna_outl[, .(strs_id, gene, region, repeat_unit,
                                n_out_rna = num(n_outliers_dbscan_global))])

rna_samp <- unique(rna_outl[, .(strs_id, os = outlier_samples_dbscan_global)])

## ---------------------------------------------------------------------------
## 2. TABELA DE STRs COM OUTLIERS
## ---------------------------------------------------------------------------
sets <- unique(rna_meta[, .(strs_id, gene, region, repeat_unit, n_out_rna)])
setorder(sets, gene, strs_id)

fwrite(sets, file.path(out_dir, "rna_outlier_sets.tsv"), sep = "\t")
cat(sprintf("STRs com outlier: %d\n", uniqueN(sets$strs_id)))

## ---------------------------------------------------------------------------
## 3. RESUMO POR GENE
## ---------------------------------------------------------------------------
genes_summary <- rna_outl[, .(
  n_strs = uniqueN(strs_id),
  n_outliers = sum(num(n_outliers_dbscan_global), na.rm = TRUE),
  gse = paste(sort(unique(gse)), collapse = ";")
), by = gene]
setorder(genes_summary, -n_strs, gene)

fwrite(genes_summary, file.path(out_dir, "rna_genes_summary.tsv"), sep = "\t")
cat(sprintf("Genes com outlier: %d\n", nrow(genes_summary)))

## ---------------------------------------------------------------------------
## 4. TABELA LONGA STR x PACIENTE (do catalogo da coorte)
## ---------------------------------------------------------------------------
cat_cols <- c("STRs_ID", "sample_id", "group", "allele1_est", "allele2_est",
              "chrom", "start", "end", "region", "repeat_unit")
catc <- fread(catalog_file, header = TRUE, sep = "\t", select = cat_cols)
cat(sprintf("Catalogo carregado: %d linhas\n", nrow(catc)))

marker <- function(ids) data.table(strs_id = unique(ids))
union_sid <- marker(rna_strs)
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

rna_os_key <- if (nrow(rna_samp)) unique(rna_samp, by = "strs_id") else data.table()

per_str_list <- lapply(unique(pat$strs_id), function(sid) {
  info <- sets[strs_id == sid, ][1]
  if (!nrow(info)) return(NULL)
  p1l <- pat[strs_id == sid]
  sc <- sum_stats(p1l[group == "case"]$maior_alelo)
  sn <- sum_stats(p1l[group == "control"]$maior_alelo)

  os_r <- if (nrow(rna_os_key)) rna_os_key[strs_id == sid]$os[1]  else NA_character_
  orr <- count_out_grp(parse_ids(os_r), p1l)

  ov <- "sem_dados"
  if (sc$n && sn$n) ov <- if (sc$min <= sn$max && sn$min <= sc$max) "sim" else "nao"

  data.table(
    strs_id = sid,
    gene = info$gene,
    region = info$region,
    repeat_unit = info$repeat_unit,
    in_rna = 1L,
    n_case = sc$n, n_control = sn$n,
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
## 6. RESUMO FINAL
## ---------------------------------------------------------------------------
cat("\n=== RESULTADO DESCRITIVO POR STR (caso x controle) ===\n")
cat(sprintf("  STRs com outlier RNA: %d\n", uniqueN(per_str$strs_id)))
cat(sprintf("  overlap maior alelo 'sim': %d | 'nao': %d | 'sem_dados': %d\n",
            sum(per_str$overlap_maior_alealo_grupos == "sim"),
            sum(per_str$overlap_maior_alealo_grupos == "nao"),
            sum(per_str$overlap_maior_alealo_grupos == "sem_dados")))
sub <- per_str[in_rna == 1L]
if (nrow(sub)) {
  wm <- function(m, n) if (sum(n) > 0) sum(m * n) / sum(n) else NaN
  cat(sprintf("  [rna] STRs=%d | overlap sim=%d nao=%d sem_dados=%d | caso n=%d (med=%.1f) control n=%d (med=%.1f)\n",
              nrow(sub),
              sum(sub$overlap_maior_alealo_grupos == "sim"),
              sum(sub$overlap_maior_alealo_grupos == "nao"),
              sum(sub$overlap_maior_alealo_grupos == "sem_dados"),
              sum(sub$n_case), wm(sub$mean_case, sub$n_case),
              sum(sub$n_control), wm(sub$mean_control, sub$n_control)))
}

cat("\n=== Saidas em:", out_dir, "===\n")
cat(paste(sort(list.files(out_dir)), collapse = "\n  "), "\n")
cat("\n=== FIM compare_gwas_rna.R ===\n")
