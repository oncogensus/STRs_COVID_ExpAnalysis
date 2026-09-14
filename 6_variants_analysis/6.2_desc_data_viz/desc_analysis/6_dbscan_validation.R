#!/usr/bin/env Rscript
# 6_dbscan_validation.R
# ---------------------------------------------------------------------------
# PROPOSITO
#   Validacao tecnica do DBSCAN: painel dual com metricas por regiao genomica.
#     Painel A: Distribuicao de genotipos (1 Cluster, 2 Clusters, 3+, Unknown)
#     Painel B: Tiers de noise (High Quality < 0.10, Acceptable < 0.25, Other)
#   Avalia TODOS os loci da coorte (sem filtro de outliers).
#
# ENTRADAS (por argumentos de linha de comando)
#   --str-catalog    STRs_analysis_dataset.tsv (coorte global)
#   --out-dir        Diretorio de saida
#
# SAIDAS
#   dbscan_dual_panel_validation.png
#   quality_funnel_summary.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ==========================================
# Parse arguments
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  return(default)
}

path_str_catalog <- parse_arg("--str-catalog")
out_dir          <- parse_arg("--out-dir", ".")

if (is.null(path_str_catalog)) {
  stop("Argumento ausente: --str-catalog")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 1. Load data
# ==========================================
cat("--- DBSCAN Technical Validation (Cohort Global) ---\n")

df_strs <- fread(path_str_catalog, header = TRUE, sep = "\t")
cat(sprintf("  STRs_analysis_dataset.tsv: %d linhas\n", nrow(df_strs)))

# ==========================================
# 2. Aggregate per STR x Sample
# ==========================================
locus_tech_stats <- df_strs[, .(
  noise_ratio = first(noise_ratio),
  n_clusters  = first(n_clusters),
  region      = first(region)
), by = .(STRs_ID, sample_id)]

# Rename regions for display
locus_tech_stats[, region := case_when(
  region == "five_prime_utr"   ~ "5' UTR",
  region == "three_prime_utr"  ~ "3' UTR",
  region == "intron"           ~ "Intron",
  region == "exon"             ~ "Exon",
  region == "promoter"         ~ "Promoter",
  region == "intergenic"       ~ "Intergenic",
  region == "non_coding_exons" ~ "Non-coding exons",
  TRUE                         ~ region
)]

# Quality tiers
locus_tech_stats[, noise_tier := case_when(
  noise_ratio < 0.10 ~ "High Quality (noise < 0.10)",
  noise_ratio < 0.25 ~ "Acceptable (noise < 0.25)",
  TRUE                ~ "Other (> 0.25)"
)]

# Cluster tiers
locus_tech_stats[, cluster_tier := case_when(
  n_clusters == 1 ~ "1 Cluster",
  n_clusters == 2 ~ "2 Clusters",
  n_clusters >= 3 ~ "3+ Clusters",
  TRUE            ~ "Unknown"
)]

# Passed QC
locus_tech_stats[, passed_qc := noise_ratio < 0.10 & n_clusters == 1]

cat(sprintf("  Observacoes unicas (STR x Sample): %d\n", nrow(locus_tech_stats)))

# ==========================================
# 3. Quality funnel summary table
# ==========================================
quality_table <- locus_tech_stats[, .(
  total_observations = .N,
  passed_qc_count    = sum(passed_qc, na.rm = TRUE),
  perc_passed        = round(sum(passed_qc, na.rm = TRUE) / .N * 100, 1),
  avg_noise          = round(mean(noise_ratio, na.rm = TRUE), 6)
), by = region]

quality_table <- quality_table[order(-total_observations)]

cat("\n--- Quality Funnel Summary ---\n")
print(as.data.frame(quality_table))

# ==========================================
# 4. Panel A: Genotypes per Sample/Region
# ==========================================
cat("\nGerando dual-panel validation plot...\n")

panel_a_data <- locus_tech_stats[, .N, by = .(region, cluster_tier)]
panel_a_data[, total := sum(N), by = region]
panel_a_data[, p := N / total]

p1 <- ggplot(panel_a_data, aes(x = region, y = p, fill = cluster_tier)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("1 Cluster"    = "#A2D2FF",
                                "2 Clusters"   = "#3498DB",
                                "3+ Clusters"  = "#2C3E50",
                                "Unknown"      = "#BDC3C7")) +
  labs(title = "A. Genotypes per Sample/Region",
       x = "", y = "% of Sample-Locus", fill = "Clusters") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 11, face = "bold"),
    plot.title  = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# ==========================================
# 5. Panel B: Noise Tiers per Sample/Region
# ==========================================
panel_b_data <- locus_tech_stats[, .N, by = .(region, noise_tier)]
panel_b_data[, total := sum(N), by = region]
panel_b_data[, p := N / total]

p2 <- ggplot(panel_b_data, aes(x = region, y = p, fill = noise_tier)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("High Quality (noise < 0.10)" = "#27AE60",
                                "Acceptable (noise < 0.25)"   = "#F1C40F",
                                "Other (> 0.25)"              = "#C0392B")) +
  labs(title = "B. Noise Tiers per Sample/Region",
       x = "", y = "% of Sample-Locus", fill = "Quality") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 11, face = "bold"),
    plot.title  = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# ==========================================
# 6. Compose and save
# ==========================================
final_layout <- p1 | p2

out_png <- file.path(out_dir, "dbscan_dual_panel_validation.png")
ggsave(
  filename = out_png,
  plot = final_layout,
  width = 14,
  height = 7,
  dpi = 600,
  bg = "white"
)
cat(sprintf("\nDual-panel salvo em: %s\n", out_png))

# Save quality table
out_csv <- file.path(out_dir, "quality_funnel_summary.csv")
write_csv2(as.data.frame(quality_table), out_csv)
cat(sprintf("Quality table salva em: %s\n", out_csv))

cat("\nConcluido.\n")
