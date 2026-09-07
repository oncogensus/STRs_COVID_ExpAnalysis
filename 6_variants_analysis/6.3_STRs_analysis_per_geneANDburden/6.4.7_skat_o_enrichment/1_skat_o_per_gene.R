# ============================================================================
# 1_skat_o_per_gene.R -- SKAT-O por gene (burden + SKAT combinados)
# ----------------------------------------------------------------------------
# Semelhante ao burden_gwas.R da 6.4.3, mas roda SKAT.O (rho otimo) em cada
# gene. Relata os 3 p-values (skat_o, skat, burden), o rho otimo e a correcao
# BH. Estrategias: gwas_burden | rna_burden | full.
#
# Uso:
#   Rscript 1_skat_o_per_gene.R [--strategy gwas_burden|rna_burden|full]
#                                [--background <arquivo>] [--out-dir <dir>]
# ============================================================================
source("00_common.R")

cmd_args <- commandArgs(trailingOnly = TRUE)
strategy       <- get_opt(cmd_args, "--strategy", "gwas_burden")
background_opt <- get_opt(cmd_args, "--background")
out_dir_opt    <- get_opt(cmd_args, "--out-dir")
min_strs       <- as.integer(get_opt(cmd_args, "--min-strs", "1"))
kernel         <- get_opt(cmd_args, "--kernel", "linear.weighted")

if (strategy == "gwas_burden") {
  out_dir <- file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_skat_o_genes")
} else if (strategy == "rna_burden") {
  out_dir <- file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_skat_o_genes_rna")
} else {
  out_dir <- file.path(REPO_ROOT,
    "6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.4.7_skat_o_enrichment/results_skat_o_genes_full")
}
if (!is.null(out_dir_opt)) out_dir <- out_dir_opt

cat("=== 1_skat_o_per_gene [strategy:", strategy, "] ===\n")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inp <- load_inputs(strategy, background_opt)
cat("Matrix:", nrow(inp$M_dosage), "amostras x", ncol(inp$M_dosage), "STRs\n")
cat("[debug] range allele2_est:", range(inp$M_dosage), "\n")
cat("[debug] genes unicos:", length(unique(inp$str_meta$gene_id)), "\n")

cat("Ajustando modelo nulo (group ~ age + sex + EV1:3)...\n")
obj <- build_null_model(inp$covar, inp$M_dosage)

genes <- unique(inp$str_meta$gene_id[!is.na(inp$str_meta$gene_id)])
cat("Testando", length(genes), "genes...\n")

res <- list()
for (g in genes) {
  strs_g <- inp$str_meta$STRs_ID[inp$str_meta$gene_id == g]
  r <- run_skat_o(strs_g, inp$M_dosage, obj,
                  min_var = min_strs, kernel = kernel)
  gname <- unique(inp$str_meta$gene_name[inp$str_meta$gene_id == g])
  res[[g]] <- data.frame(
    gene = g,
    gene_name = ifelse(length(gname) > 0, gname[1], NA),
    strs_ids = paste(strs_g, collapse = ";"),
    r, stringsAsFactors = FALSE)
}

skat <- bind_rows(res)
skat$q_value_skat_o   <- p.adjust(skat$p_value_skat_o, method = "BH")
skat$q_value_skat     <- p.adjust(skat$p_value_skat, method = "BH")
skat$q_value_burden   <- p.adjust(skat$p_value_burden, method = "BH")
skat <- skat[order(skat$p_value_skat_o), ]
if (anyNA(skat$p_value_skat_o)) skat <- skat[order(!is.na(skat$p_value_skat_o),
                                                   skat$p_value_skat_o), ]

write_hits(skat, "p_value_skat_o", "q_value_skat_o",
           "skat_o_per_gene", out_dir)
write_hits(skat, "p_value_burden", "q_value_burden",
           "burden_per_gene", out_dir)

cat("Resultados em:", out_dir, "\n")
cat("\n=== SKAT-O por gene (top 20) ===\n")
print(head(skat, 20))
unc <- skat[!is.na(skat$q_value_skat_o) & skat$q_value_skat_o < 0.05, ]
cat("\n=== SKAT-O hits corrigidos (q<0.05):", nrow(unc), "===\n")
print(unc)
cat("\n=== FIM 1_skat_o_per_gene.R ===\n")