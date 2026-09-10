#!/usr/bin/env Rscript
# 6.2.4_str_coverage_per_patient.R
# ---------------------------------------------------------------------------
# PURPOSE
#   Computes, PER PATIENT, the genomic coverage of STRs called after STRling
#   quality control (the QC that generates STRs_analysis_dataset).
#   For each sample:
#     n_strs        = number of distinct STR loci
#     bp_covered    = covered bases (overlapping intervals merged per chr,
#                     no double counting)
#     genome_frac   = bp_covered / genome size (chr 1-22, X, Y)
#   Produces descriptive statistics (mean, max, median, min, SD, % genome)
#   globally and per group (case/control).
#
# INPUTS
#   --dataset   STRs_analysis_dataset.tsv  (already STRling-QC filtered)
#   --genome    genome.txt (chr<TAB>size) [default: 6.2_desc_data_viz/desc_analysis/genome.txt]
#   --out-dir   Output directory (created if missing)
#   --valid-genotype-only  If set, keep only STRs with allele1_est AND
#                          allele2_est non-NA (otherwise uses everything post-QC).
#
# OUTPUTS
#   coverage_by_patient.tsv   (one row per sample)
#   coverage_summary.tsv      (descriptive stats, global + per group)
#
# Usage:
#   Rscript 6.2.4_str_coverage_per_patient.R \
#     --dataset $REPO/samples/STRs_analysis_dataset.tsv \
#     --genome  $REPO/6_variants_analysis/6.2_desc_data_viz/desc_analysis/genome.txt \
#     --out-dir results/coverage
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table) })

get_opt <- function(args, flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
has_flag <- function(args, flag) flag %in% args

args <- commandArgs(trailingOnly = TRUE)
dataset_file <- get_opt(args, "--dataset", NULL)
genome_file  <- get_opt(args, "--genome", NULL)
out_dir      <- get_opt(args, "--out-dir", "results_coverage")

if (is.null(dataset_file) || is.null(genome_file))
  stop("Usage: Rscript 6.2.4_str_coverage_per_patient.R --dataset <tsv> --genome <genome.txt> [--out-dir <dir>] [--valid-genotype-only]")
for (f in c(dataset_file, genome_file)) if (!file.exists(f)) stop("missing file: ", f)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== STR genomic coverage per patient (post-QC) ===\n")
cat("dataset   :", dataset_file, "\n")
cat("genome    :", genome_file, "\n")
cat("out_dir   :", out_dir, "\n")

## ---------------------------------------------------------------
## 1. GENOME (chr 1-22, X, Y) -> total size in bp
## ---------------------------------------------------------------
gen <- fread(genome_file, header = FALSE, sep = "\t")
setnames(gen, c("chrom", "size"))
gen <- gen[grepl("^chr([0-9]{1,2}|X|Y)$", chrom)]
gen[, chrom := sub("^chr", "", chrom)]
genome_total <- sum(gen$size)
cat(sprintf("Total genome (chr 1-22, X, Y): %s bp (%.1f Mb)\n",
            formatC(genome_total, big.mark = ","), genome_total / 1e6))

## ---------------------------------------------------------------
## 2. DATASET (already STRling-QC filtered)
## ---------------------------------------------------------------
strs <- fread(dataset_file, header = TRUE, sep = "\t")
nec <- c("sample_id", "group", "chrom", "start", "end")
if (!all(nec %in% names(strs)))
  stop("dataset must contain: ", paste(nec, collapse = ", "))
strs[, sample_id := as.character(sample_id)]
strs[, group := tolower(trimws(as.character(group)))]
strs[, chrom := sub("^chr", "", as.character(chrom))]
strs[, start := as.numeric(start)][, end := as.numeric(end)]
strs <- strs[!is.na(start) & !is.na(end)]


n_pat_total <- uniqueN(strs$sample_id)
cat(sprintf("Samples in dataset: %d | STR x patient rows: %s\n",
            n_pat_total, formatC(nrow(strs), big.mark = ",")))

## ---------------------------------------------------------------
## 3. COVERAGE PER PATIENT (intervals merged per chr)
## ---------------------------------------------------------------
# merge intervals within a chromosome -> total length (bp)
merge_chr_width <- function(dd) {
  if (!nrow(dd)) return(0L)
  dd <- dd[order(start, end)]
  ss <- dd$start; ee <- dd$end
  is <- ss[1]; ie <- ee[1]; tot <- 0L
  if (length(ss) > 1) {
    for (i in seq.int(2, length(ss))) {
      if (ss[i] <= ie) {
        if (ee[i] > ie) ie <- ee[i]
      } else {
        tot <- tot + (ie - is + 1L)
        is <- ss[i]; ie <- ee[i]
      }
    }
  }
  tot + (ie - is + 1L)
}

samples <- unique(strs$sample_id)
cols_out <- data.table(
  sample_id = character(), group = character(),
  n_strs = integer(), bp_covered = numeric(),
  genome_frac = numeric()
)

for (sid in samples) {
  sub <- strs[sample_id == sid, .(chrom, start, end)]
  sub <- unique(sub)                          # distinct loci (chrom,start,end)
  n_strs <- nrow(sub)
  bp <- 0
  if (n_strs && !all(is.na(sub$start))) {
    bp <- sum(vapply(split(sub, by = "chrom", keep.by = FALSE),
                     merge_chr_width, numeric(1)))
  }
  cols_out <- rbind(cols_out, data.table(
    sample_id = sid,
    group = strs[sample_id == sid, group[1]],
    n_strs = as.integer(n_strs),
    bp_covered = as.numeric(bp),
    genome_frac = bp / genome_total
  ))
}
cat(sprintf("Coverage computed for %d patients.\n", nrow(cols_out)))

## ---------------------------------------------------------------
## 4. DESCRIPTIVE STATISTICS
## ---------------------------------------------------------------
desc_one <- function(dt) {
  data.table(
    n_patients = nrow(dt),
    n_strs_mean = mean(dt$n_strs), n_strs_sd = sd(dt$n_strs),
    n_strs_median = median(dt$n_strs), n_strs_min = min(dt$n_strs),
    n_strs_max = max(dt$n_strs),
    bp_mean = mean(dt$bp_covered), bp_sd = sd(dt$bp_covered),
    bp_max = max(dt$bp_covered),
    bp_median = median(dt$bp_covered), bp_min = min(dt$bp_covered),
    genome_frac_mean = mean(dt$genome_frac), genome_frac_sd = sd(dt$genome_frac),
    genome_frac_max = max(dt$genome_frac),
    genome_frac_median = median(dt$genome_frac),
    genome_frac_min = min(dt$genome_frac),
    genome_pct_mean = mean(dt$genome_frac) * 100
  )
}

d_global <- desc_one(cols_out)
d_global[, group := "GLOBAL"]
d_group <- cols_out[, desc_one(.SD), by = group]

summary_tbl <- rbindlist(list(d_group, d_global), use.names = TRUE, fill = TRUE)
setcolorder(summary_tbl, c("group", setdiff(names(summary_tbl), "group")))

## ---------------------------------------------------------------
## 5. OUTPUTS
## ---------------------------------------------------------------
fwrite(cols_out, file.path(out_dir, "coverage_by_patient.tsv"), sep = "\t")
fwrite(summary_tbl, file.path(out_dir, "coverage_summary.tsv"), sep = "\t")

## ---------------------------------------------------------------
## 6. STDOUT SUMMARY
## ---------------------------------------------------------------
nice_pct <- function(x, digs = 3) sprintf(paste0("%.", digs, "f%%"), x)
cat("\n=== DESCRIPTIVE (per patient, post-QC) ===\n")
cat(sprintf("Global: n=%d | n_strs mean=%.0f max=%d | bp mean=%.0f max=%.0f | %%genome mean=%s max=%s\n",
            d_global$n_patients, d_global$n_strs_mean, d_global$n_strs_max,
            d_global$bp_mean, d_global$bp_max,
            nice_pct(d_global$genome_pct_mean),
            nice_pct(d_global$genome_frac_max * 100)))
if (nrow(d_group)) for (i in seq_len(nrow(d_group))) {
  g <- d_group[i]
  cat(sprintf("Group %s: n=%d | n_strs mean=%.0f max=%d | bp mean=%.0f max=%.0f | %%genome mean=%s\n",
              g$group, g$n_patients, g$n_strs_mean, g$n_strs_max,
              g$bp_mean, g$bp_max, nice_pct(g$genome_pct_mean)))
}
cat("\n=== Outputs in:", out_dir, "===\n")
cat(paste(sort(list.files(out_dir)), collapse = "\n  "), "\n")
cat("\n=== END 6.2.4_str_coverage_per_patient.R ===\n")