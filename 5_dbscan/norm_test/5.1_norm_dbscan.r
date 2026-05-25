### STR Normalization using Linear Regression for DBSCAN Analysis
library(data.table)
library(dplyr)
library(stringr)

# ==========================================
# 0. Global Normalization Function
# ==========================================
# This ensures consistency across all files
normalize_ids <- function(ids) {
  ids %>%
    as.character() %>%
    toupper() %>%
    # 1. Remove .bam and anything following it
    str_remove("(?i)[._-]?\\d*BAM.*$") %>%
    # 2. FIXED: Remove replicate suffix (-1) only if preceded by a digit
    # Prevents GIG-MT-4 from becoming GIG-MT
    str_remove("(?<=\\d)-[0-9]$") %>%
    # 3. Standardize leading zeros (REC-07 -> REC-7)
    str_replace_all("-0+([0-9]+)", "-\\1") %>%
    str_trim()
}

# ==========================================
# 1. Load STR data
# ==========================================
cat("\nSTEP 1: Reading STRs...\n")
path_strs <- "../../samples/STRs_annotated_region.tsv"
dat.strs <- fread(path_strs, header=T, sep="\t", data.table=F)

dat.strs$sample_id_clean <- normalize_ids(dat.strs$sample_id)

cat("Sample ID example (Original):", head(unique(dat.strs$sample_id), 2), "\n")
cat("Sample ID example (Cleaned): ", head(unique(dat.strs$sample_id_clean), 2), "\n")

# ==========================================
# 2. Load PCAs (Ancestry Coordinates)
# ==========================================
cat("\nSTEP 2: Reading PCA coordinates...\n")
path_pca <- "../../4_ancestry/EthSEQ_Results_3D/Report.PCAcoord"
dat.pca <- read.delim(path_pca, header=T, sep="\t", stringsAsFactors=F)

# Standardize column name and clean IDs
names(dat.pca)[names(dat.pca) == "sample.id"] <- "sample_id" 
dat.pca$sample_id_clean <- normalize_ids(dat.pca$sample_id)

cat("PCA ID example (Cleaned): ", head(unique(dat.pca$sample_id_clean), 2), "\n")

# ==========================================
# 3. Load Demographic Info (Phenotypes)
# ==========================================
cat("\nSTEP 3: Reading phenotypes...\n")
path_pheno <- "../../samples/samples_infos.csv"
dat.pheno <- fread(path_pheno, header=T, sep=",", data.table=F)

# Detect sample column and clean
# (Adjust 'sample' to the actual column name in your CSV if different)
dat.pheno$sample_id_clean <- normalize_ids(dat.pheno$sample)

dat.pheno.final <- dat.pheno %>%
  select(sample_id_clean, age, sex, group) %>%
  mutate(
    age = as.numeric(age),
    sex = as.factor(sex)
  )

cat("Pheno ID example (Cleaned):", head(unique(dat.pheno.final$sample_id_clean), 2), "\n")

# ==========================================
# 4. Merge Data for Regression
# ==========================================
cat("\nSTEP 4: Merging all data sources...\n")

ids_strs  <- unique(dat.strs$sample_id_clean)
ids_pca   <- unique(dat.pca$sample_id_clean)
ids_pheno <- unique(dat.pheno.final$sample_id_clean)

common_ids <- intersect(intersect(ids_strs, ids_pca), ids_pheno)
cat("Samples common to all files:", length(common_ids), "/ 168 expected\n")

if(length(common_ids) < length(ids_strs)) {
    cat("⚠️ Warning: Some STR samples are missing PCA or Pheno data.\n")
}

# Final merged dataset for regression
dat.merged <- dat.strs %>%
  inner_join(dat.pca %>% select(sample_id_clean, EV1, EV2, EV3), by="sample_id_clean") %>%
  inner_join(dat.pheno.final, by="sample_id_clean")

# Integrity filter
dat.final <- dat.merged %>%
  filter(
    !is.na(allele2_est), !is.na(depth),
    !is.na(EV1), !is.na(EV2), !is.na(EV3),
    !is.na(sex), !is.na(age)
  )

cat("Final row count for regression:", nrow(dat.final), "\n")

# ==========================================
# 5. Normalization by Linear Regression
# ==========================================
cat("\nSTEP 5: Running regression-based normalization...\n")

if(nrow(dat.final) == 0) stop("FATAL: Dataframe is empty. Check your merges/filters.")

dat.normalized <- dat.final %>%
  group_by(STRs_ID) %>%
  mutate(
    n_samples = n(), 
    # Only run regression if we have enough data points (n > 10)
    allele2_residuals = if(n() > 10) { 
      tryCatch({
        mod <- lm(allele2_est ~ sex + EV1 + EV2 + EV3 + depth + age, data=cur_data())
        residuals(mod)
      }, error = function(e) NA_real_)
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

# ==========================================
# 6. Save Results
# ==========================================
cat("\nSTEP 6: Saving normalized residuals...\n")

dat.out <- dat.normalized %>%
  select(
    STRs_ID, sample_id, sample_id_clean, group, 
    region, gene_id, gene_name,
    depth, allele2_est, allele2_residuals
  )

fwrite(dat.out, "STRs_normalized_residuals.tsv", sep="\t")
cat("✔ Success! File saved: STRs_normalized_residuals.tsv\n")

# ==========================================
# 7. Summary Statistics
# ==========================================
cat("\n=== FINAL STATISTICS ===\n")
cat("Unique STRs processed: ", length(unique(dat.normalized$STRs_ID)), "\n")
cat("Unique Patients processed: ", length(unique(dat.normalized$sample_id_clean)), "\n")

res_valid <- dat.normalized$allele2_residuals[!is.na(dat.normalized$allele2_residuals)]
if(length(res_valid) > 0) {
    cat("Normalization Summary (Residuals):\n")
    print(summary(res_valid))
}