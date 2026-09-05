# compare_burden_hits.R
# ---------------------------------------------------------------------------
# Compara os hits do burden entre as estrategias gwas_burden e rna_burden:
#   - genes em sobreposicao (uncorrected p<0.05 e corrected BH q<0.05)
#   - uniao completa com flags hit_gwas / hit_rna e os p/q de cada lado
#
# Entradas (padroes cluster):
#   <repo>/6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.3_burden_test/
#       results_gwas_burden/skat_per_gene_hits_{uncorrected,corrected}.tsv
#       results_rna/skat_per_gene_hits_{uncorrected,corrected}.tsv
#
# Saidas (--out-dir, padrao results_burden_comparison/):
#   burden_hits_overlap_uncorrected.tsv / burden_hits_overlap_corrected.tsv
#   burden_hits_union_uncorrected.tsv   / burden_hits_union_corrected.tsv
#
# Uso:
#   Rscript compare_burden_hits.R [--repo <repodir>] [--out-dir <dir>]
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table) })

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

cmd_args <- commandArgs(trailingOnly = TRUE)
REPO_ROOT <- get_opt(cmd_args, "--repo",
                     "/storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis")

burden_dir <- file.path(REPO_ROOT,
  "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.3_burden_test")
out_dir <- get_opt(cmd_args, "--out-dir",
                   file.path(burden_dir, "results_burden_comparison"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hit_files <- list(
  gwas_unc = file.path(burden_dir, "results_gwas_burden/skat_per_gene_hits_uncorrected.tsv"),
  gwas_cor = file.path(burden_dir, "results_gwas_burden/skat_per_gene_hits_corrected.tsv"),
  rna_unc  = file.path(burden_dir, "results_rna/skat_per_gene_hits_uncorrected.tsv"),
  rna_cor  = file.path(burden_dir, "results_rna/skat_per_gene_hits_corrected.tsv")
)

read_hits <- function(path) {
  if (!file.exists(path)) stop("arquivo ausente: ", path)
  dt <- fread(path, header = TRUE, sep = "\t", data.table = FALSE)
  if (!"gene" %in% colnames(dt)) stop("coluna 'gene' ausente em: ", path)
  dt[!is.na(dt$p_value), , drop = FALSE]
}

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

  # Compartilhados primeiro, depois por p (RNA) e por p (GWAS); NAs vao para o fim
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

g_unc <- read_hits(hit_files$gwas_unc)
g_cor <- read_hits(hit_files$gwas_cor)
r_unc <- read_hits(hit_files$rna_unc)
r_cor <- read_hits(hit_files$rna_cor)

add_flag <- function(dt, name) {
  dt[[name]] <- rep(1L, nrow(dt))
  dt
}
g_unc <- add_flag(g_unc, "hit_gwas")
g_cor <- add_flag(g_cor, "hit_gwas")
r_unc <- add_flag(r_unc, "hit_rna")
r_cor <- add_flag(r_cor, "hit_rna")

u_unc <- build_union(g_unc, r_unc, "hit_gwas", "hit_rna")
u_cor <- build_union(g_cor, r_cor, "hit_gwas", "hit_rna")

ov_unc <- u_unc[u_unc$hit_gwas == 1L & u_unc$hit_rna == 1L, , drop = FALSE]
ov_cor <- u_cor[u_cor$hit_gwas == 1L & u_cor$hit_rna == 1L, , drop = FALSE]

fwrite(u_unc, file.path(out_dir, "burden_hits_union_uncorrected.tsv"), sep = "\t")
fwrite(u_cor, file.path(out_dir, "burden_hits_union_corrected.tsv"), sep = "\t")
fwrite(ov_unc, file.path(out_dir, "burden_hits_overlap_uncorrected.tsv"), sep = "\t")
fwrite(ov_cor, file.path(out_dir, "burden_hits_overlap_corrected.tsv"), sep = "\t")

cat("=== Comparativo burden GWAS x RNA ===\n")
cat("  uncorrected (p<0.05): GWAS=", nrow(g_unc), " RNA=", nrow(r_unc),
    " overlap=", nrow(ov_unc), "\n", sep = "")
cat("  corrected   (q<0.05): GWAS=", nrow(g_cor), " RNA=", nrow(r_cor),
    " overlap=", nrow(ov_cor), "\n", sep = "")
cat("Saidas em:", out_dir, "\n")

if (nrow(ov_unc) > 0) {
  cat("\n=== Genes em sobreposicao (uncorrected) ===\n")
  print(ov_unc[, c("gene", "gene_name", "gwas_n_variants", "rna_n_variants",
                   "gwas_p_value", "rna_p_value")])
}