cat("\n========================================\n")
cat("ANCESTRY x OUTLIERS HIGH RESOLUTION\n")
cat("========================================\n")

# Loading essential libraries
library(data.table)
library(rstatix)
library(stringr)

# =========================
# 0. PATHS
# =========================

input_file <- "../../samples/STRs_analysis_dataset.tsv"
output_dir <- "results/high_resolution"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("[DEBUG] Created directory:", output_dir, "\n")
}

# =========================
# 1. LOAD DATA
# =========================

cat("\nLoading data...\n")
merged_dt <- fread(input_file, sep = "\t")

cat("[DEBUG] Initial Dataset Dimensions:", nrow(merged_dt), "rows,", ncol(merged_dt), "columns\n")

# =========================
# 2. QUALITY CONTROL (QC) CRITERIA
# =========================

cat("\nApplying DBSCAN Quality Control...\n")

# NEW CRITERIA: A record passes QC if it has clusters, at least one outlier, 
# and the noise ratio is 10% or less.
merged_dt[, dbscan_pass := (
  n_clusters > 0 &
  n_outliers >= 1 &
  noise_ratio <= 0.10
)]

cat("[DEBUG] QC Pass Rate:", round(mean(merged_dt$dbscan_pass) * 100, 2), "%\n")

# =========================
# 3. BUILD ANCESTRY (GENERALIZED)
# =========================

cat("\nBuilding ancestry proportions...\n")

all_ancestries <- unique(
  unlist(
    str_extract_all(
      merged_dt$contribution,
      "[A-Z]+(?=\\()"
    )
  )
)
all_ancestries <- all_ancestries[!is.na(all_ancestries)]

parse_contribution <- function(contrib_string, ancestries) {
  res <- setNames(rep(0, length(ancestries)), ancestries)
  if (is.na(contrib_string) || contrib_string == "") return(res)
  parts <- strsplit(contrib_string, "\\|")[[1]]
  for (p in parts) {
    p <- trimws(p)
    match <- regmatches(p, regexec("([A-Z]+)\\(([^%)]+)%?\\)", p))[[1]]
    if (length(match) == 3) {
      pop <- match[2]
      value <- suppressWarnings(as.numeric(match[3])) / 100
      if (!is.na(value) && pop %in% names(res)) res[pop] <- value
    }
  }
  return(res)
}

ancestry_dt <- unique(merged_dt[, .(
  sample_id,
  type = type,
  pop = pop,
  contribution = contribution
)])

for (anc in all_ancestries) ancestry_dt[, (anc) := 0]
for (anc in all_ancestries) ancestry_dt[type == "INSIDE" & pop == anc, (anc) := 1]

closest_idx <- which(ancestry_dt$type == "CLOSEST")
if (length(closest_idx) > 0) {
  parsed_list <- lapply(ancestry_dt$contribution[closest_idx], function(x) as.list(parse_contribution(x, all_ancestries)))
  parsed <- rbindlist(parsed_list, fill = TRUE)
  ancestry_dt[closest_idx, (all_ancestries) := parsed]
  
  ancestry_dt[closest_idx, total := Reduce(`+`, .SD), .SDcols = all_ancestries]
  ancestry_dt[closest_idx, (all_ancestries) := lapply(.SD, function(x) x / total), .SDcols = all_ancestries]
}

# =========================
# 4. COMBINE ANCESTRY
# =========================

cat("\nConsolidating one row per sample...\n")
samples_with_inside <- unique(ancestry_dt[type == "INSIDE", sample_id])

final_ancestry <- data.table(sample_id = character())
for (anc in all_ancestries) set(final_ancestry, j = anc, value = numeric())

if (length(samples_with_inside) > 0) {
  inside_part <- ancestry_dt[type == "INSIDE" & sample_id %in% samples_with_inside, .SD[1], by = sample_id]
  final_ancestry <- rbind(final_ancestry, inside_part[, c("sample_id", all_ancestries), with = FALSE])
}

samples_without_inside <- setdiff(unique(ancestry_dt$sample_id), samples_with_inside)
if (length(samples_without_inside) > 0) {
  closest_part <- ancestry_dt[type == "CLOSEST" & sample_id %in% samples_without_inside, .SD[1], by = sample_id]
  final_ancestry <- rbind(final_ancestry, closest_part[, c("sample_id", all_ancestries), with = FALSE])
}

# =========================
# 5. MERGE BACK
# =========================

cat("\nMerging ancestry back to main dataset...\n")
cols_to_remove <- intersect(colnames(merged_dt), all_ancestries)
if(length(cols_to_remove) > 0) merged_dt[, (cols_to_remove) := NULL]
merged_dt <- merge(merged_dt, final_ancestry, by = "sample_id", all.x = TRUE)

# =========================
# 6. OUTLIERS & STRENGTH (Residuals)
# =========================

cat("\nFlagging outliers and extracting residuals (strength)...\n")

merged_dt[, is_outlier := str_detect(outlier_samples, fixed(sample_id))]
merged_dt[is.na(is_outlier) | outlier_samples == "", is_outlier := FALSE]

merged_dt[, outlier_strength_val := as.numeric(outlier_residuals)]
merged_dt[is_outlier == FALSE, outlier_strength_val := 0]
merged_dt[is_outlier == TRUE & is.na(outlier_strength_val), outlier_strength_val := 0]

# =========================
# 7. AGGREGATE
# =========================

cat("\nAggregating by Region and Sample (QC Filtered - Main Dataset)...\n")

# LÓGICA DO BK: Utilizada para gerar o dataset principal rigoroso
region_sample_dt <- merged_dt[, .(
  total_STRs_raw = .N,
  total_STRs_after_QC = sum(dbscan_pass),
  n_outliers = sum(is_outlier & dbscan_pass),
  outlier_prop = sum(is_outlier & dbscan_pass) / max(1, sum(dbscan_pass)),
  outlier_strength = sum(abs(outlier_strength_val[dbscan_pass]), na.rm = TRUE) / max(1, sum(is_outlier & dbscan_pass))
), by = c("sample_id", "region")]

region_sample_dt[is.nan(outlier_strength), outlier_strength := 0]
region_sample_dt <- merge(region_sample_dt, final_ancestry, by = "sample_id", all.x = TRUE)


cat("\nAggregating by Region and Sample (Unfiltered - ONLY for Correlation/Debug)...\n")

# LÓGICA NOVA: Objeto à parte sem filtro DBSCAN para alimentar apenas a correlação e o Debug
region_sample_dt_unfiltered <- merged_dt[, .(
  total_STRs_raw = .N,
  total_STRs_after_QC = sum(dbscan_pass, na.rm = TRUE),
  n_outliers = sum(is_outlier, na.rm = TRUE),
  outlier_prop = mean(is_outlier, na.rm = TRUE),
  outlier_strength = sum(abs(outlier_strength_val), na.rm = TRUE) / max(1, sum(is_outlier, na.rm = TRUE))
), by = c("sample_id", "region")]

region_sample_dt_unfiltered[is.nan(outlier_strength), outlier_strength := 0]
region_sample_dt_unfiltered <- merge(region_sample_dt_unfiltered, final_ancestry, by = "sample_id", all.x = TRUE)

# ==========================================
# DEBUG: INVESTIGATING MISSING REGIONS (USING UNFILTERED DATA)
# ==========================================

cat("\n--- DEBUG START (Based on Unfiltered Object) ---\n")

raw_regions <- unique(merged_dt$region)
cat("[1] Regions found in raw input:", paste(raw_regions, collapse=", "), "\n")

agg_regions_unfiltered <- unique(region_sample_dt_unfiltered$region)
cat("[2] Regions present after unfiltered aggregation:", paste(agg_regions_unfiltered, collapse=", "), "\n")

missing_targets <- c("CDS", "five_prime_utr")
cat("[3] Checking variance for target regions:\n")

for(target in missing_targets) {
  if(target %in% agg_regions_unfiltered) {
    sub_dt <- region_sample_dt_unfiltered[region == target]
    cat("  -> Region:", target, "\n")
    cat("     - Variance (outlier_prop):", var(sub_dt$outlier_prop, na.rm=TRUE), "\n")
    cat("     - Variance (outlier_strength):", var(sub_dt$outlier_strength, na.rm=TRUE), "\n")
    cat("     - Total Rows:", nrow(sub_dt), "\n")
    cat("     - Samples with outliers:", sum(sub_dt$n_outliers > 0, na.rm=TRUE), "\n")
  } else {
    cat("  -> Warning:", target, "NOT FOUND in unfiltered aggregated data.\n")
  }
}

cat("[4] Exact names of detected regions (Top 10):", paste(head(agg_regions_unfiltered, 10), collapse=", "), "\n")
cat("--- DEBUG END ---\n")

# =========================
# 8. CORRELATION
# =========================

cat("\nRunning Spearman correlations for Prop and Strength (Using Unfiltered Data)...\n")

target_metrics <- c("outlier_prop", "outlier_strength")

# Using the unrestricted object created specifically for corr_data
corr_data <- region_sample_dt_unfiltered[total_STRs_after_QC > 0]

cor_region <- corr_data[, {
  res <- list()
  for (m in target_metrics) {
    for (anc in all_ancestries) {
      x <- .SD[[anc]]
      y <- .SD[[m]]
      
      # Verificação de variância segura da versão nova
      if (length(x) > 2 && var(x, na.rm=TRUE) > 0 && var(y, na.rm=TRUE) > 0) {
        test <- suppressWarnings(cor.test(y, x, method = "spearman"))
        res[[paste0(m, "_", anc)]] <- data.table(
          estimate = test$estimate,
          p = test$p.value,
          metric = m,
          ancestry = anc
        )
      }
    }
  }
  rbindlist(res)
}, by = region]

if(nrow(cor_region) > 0) {
  cor_region[, p_adj := p.adjust(p, method = "BH")]
  cat("[DEBUG] Correlations computed successfully.\n")
} else {
  cat("[WARNING] No correlations computed. Check if ancestry columns vary.\n")
}

# =========================
# 9. SAVE ALL FILES
# =========================

cat("\nSaving results...\n")

anc_region <- merged_dt[, lapply(.SD, mean, na.rm = TRUE), by = region, .SDcols = all_ancestries]

# Export data
fwrite(region_sample_dt, file.path(output_dir, "plotdata_region_sample.csv"))
fwrite(cor_region, file.path(output_dir, "correlation_full.csv"))
fwrite(anc_region, file.path(output_dir, "ancestry_region_distribution_wide.csv"))

cat("\n========================================\n")
cat("PIPELINE COMPLETED SUCCESSFULLY\n")
cat("Files saved in:", output_dir, "\n")
cat("========================================\n")