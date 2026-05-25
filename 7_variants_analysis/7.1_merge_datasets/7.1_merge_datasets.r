library(data.table)
library(dplyr)
library(stringr)

# ==========================================
# STEP 1: Loading Datasets
# ==========================================
cat("\n", strrep("=", 42), "\n")
cat("STEP 1: Loading Datasets\n")
cat(strrep("=", 42), "\n")

# Paths based on project structure
path_annot  <- "../../samples/STRs_annotated_region.tsv"
path_dbscan <- "../../5_dbscan/outliers_search/results_dbscan/outliers_per_str.tsv"
path_eth    <- "../../4_ancestry/EthSEQ_Results_3D/Report.txt"
path_groups <- "../../samples/samples_infos.csv"

# Check if all files exist
if(!all(file.exists(c(path_annot, path_dbscan, path_eth, path_groups)))) {
  stop("One or more input files were not found. Verify your relative paths.")
}

dat.annot  <- fread(path_annot, header = TRUE)
dat.dbscan <- fread(path_dbscan, header = TRUE)
dat.eth    <- fread(path_eth, header = TRUE, sep = "auto")
dat.groups <- fread(path_groups, header = TRUE)

cat("✔ All files loaded successfully.\n")

# ==========================================
# STEP 2: Standardizing Patient IDs
# ==========================================
cat("\n", strrep("=", 42), "\n")
cat("STEP 2: Standardizing IDs\n")
cat(strrep("=", 42), "\n")

normalize_ids <- function(ids) {
  ids %>%
    as.character() %>%
    toupper() %>%
    str_remove("(?i)[._-]?\\d*BAM.*$") %>%
    str_remove("(?<=\\d)-[0-9]$") %>%
    str_replace_all("-0+([0-9]+)", "-\\1") %>%
    str_trim()
}

# Standardize column names and IDs across all frames
eth_id_col <- grep("sample|id", colnames(dat.eth), ignore.case = TRUE, value = TRUE)[1]
setnames(dat.eth, eth_id_col, "sample_id_clean")

dat.eth[, sample_id_clean := normalize_ids(sample_id_clean)]
dat.annot[, sample_id_clean := normalize_ids(sample_id)]
dat.groups[, sample_id_clean := normalize_ids(sample)]

# ==========================================
# STEP 3: Debugging Data Coverage & Counts
# ==========================================
cat("\n", strrep("=", 42), "\n")
cat("STEP 3: Debugging Data Coverage\n")
cat(strrep("=", 42), "\n")

# Extract unique identifiers for comparison
ids_annot  <- unique(dat.annot$sample_id_clean)
ids_groups <- unique(dat.groups$sample_id_clean)
ids_eth    <- unique(dat.eth$sample_id_clean)
strs_annot  <- unique(dat.annot$STRs_ID)
strs_dbscan <- unique(dat.dbscan$STRs_ID)

# Sample Level Report
cat("--- Sample Coverage (Patient Level) ---\n")
cat(sprintf("Total Patients in Annotation:    %d\n", length(ids_annot)))
cat(sprintf("Total Patients in Group Mapping: %d\n", length(ids_groups)))
cat(sprintf("Total Patients in EthSEQ:        %d\n", length(ids_eth)))

match_groups <- intersect(ids_annot, ids_groups)
match_eth    <- intersect(ids_annot, ids_eth)

cat(sprintf("\n✔ Match Annotation x Groups: %d/%d (%.1f%%)\n", 
            length(match_groups), length(ids_annot), (length(match_groups)/length(ids_annot))*100))
cat(sprintf("✔ Match Annotation x EthSEQ: %d/%d (%.1f%%)\n", 
            length(match_eth), length(ids_annot), (length(match_eth)/length(ids_annot))*100))

# Variant Level Report
cat("\n--- Variant Coverage (Locus Level) ---\n")
cat(sprintf("Total STRs in Annotation: %d\n", length(strs_annot)))
cat(sprintf("Total STRs in DBSCAN:     %d\n", length(strs_dbscan)))

match_strs <- intersect(strs_annot, strs_dbscan)
cat(sprintf("✔ Match Annotation x DBSCAN: %d/%d (%.1f%%)\n", 
            length(match_strs), length(strs_annot), (length(match_strs)/length(strs_annot))*100))

# Check for critical missing data
if(length(match_groups) < length(ids_annot)) {
  missing <- setdiff(ids_annot, ids_groups)
  warning(sprintf("%d samples are missing group metadata (Example: %s)", 
                  length(missing), paste(head(missing, 3), collapse=", ")))
}

# ==========================================
# STEP 4: Integrating Data (Merging)
# ==========================================
cat("\n", strrep("=", 42), "\n")
cat("STEP 4: Integrating Data (Merging)\n")
cat(strrep("=", 42), "\n")

dat.final <- dat.annot %>%
  left_join(dat.groups %>% select(sample_id_clean, group, age, sex), by = "sample_id_clean") %>%
  left_join(dat.eth %>% select(sample_id_clean, pop, contribution, type), by = "sample_id_clean") %>%
  left_join(dat.dbscan %>% select(STRs_ID, n_outliers, outlier_samples, outlier_residuals, n_clusters, noise_ratio), 
            by = "STRs_ID")

# ==========================================
# STEP 5: Final Export
# ==========================================
output_file <- "../../samples/STRs_analysis_dataset.tsv"

dat.final <- dat.final %>%
  select(
    # IDs & Clinical Info
    STRs_ID, group, age, sex,
    
    # Quantitative STR Metrics
    allele1_est, allele2_est, depth,
    
    # Genomic Context
    repeat_unit, gene_id, gene_name, region, chrom, start, end, sample_id,
    
    # Population Genetics
    pop, contribution, type,
    
    # Outlier Metrics
    n_outliers, outlier_samples, outlier_residuals, n_clusters, noise_ratio
  )

fwrite(dat.final, output_file, sep="\t")

cat("\n", strrep("=", 42), "\n")
cat("PROCESS COMPLETE!\n")
cat("Final dataset generated with", nrow(dat.final), "rows.\n")
cat("Output path: ", output_file, "\n")
cat(strrep("=", 42), "\n")