library(readr)
library(dplyr)
library(stringr)


# ============================================================================
# PATH CONFIGURATION 
# ============================================================================


# Main input file path
path <- "../samples/gp_global/merged/STRs_annotated_complete.tsv"


# Output file path
output_path <- "../samples/gp_global/merged/merged_summary_dataViz.tsv"

# ============================================================================
# PROCESSING
# ============================================================================


# Read entire file text
lines <- readLines(path)


# Modify only the first line to remove the '#'
lines[1] <- sub("^#", "", lines[1])


# Write to temporary file
temp_file <- tempfile(fileext = ".tsv")
writeLines(lines, temp_file)


# Read temporary file with correct header
df <- read_tsv(temp_file, col_names = TRUE)


# Proceed with extraction and calculation
df <- df %>%
  mutate(STRLEN = as.numeric(str_extract(INFO, "(?<=STRLEN=)\\d+")),
         RU = str_extract(INFO, "(?<=RU=)[^;]+"),
         RU_length = nchar(RU),
         STRLEN_x_RU_len = STRLEN * RU_length)


# Calculate summary - one line per ID
# Including Region_Type and Gene_Name (taking the first value from each group)
summary_df <- df %>%
  group_by(ID) %>%
  summarise(mean_STRLEN_x_RU_len = mean(STRLEN_x_RU_len, na.rm=TRUE),
            p_adj_flag = as.integer(any(p_adj < 0.05, na.rm=TRUE)),
            Region_Type = first(Region_Type),
            Gene_Name = first(Gene_Name),
            .groups = 'drop')


# Save output
write_tsv(summary_df, output_path)


cat("Summary file saved at:", output_path, "\n")
cat("Total rows:", nrow(summary_df), "\n")
cat("Columns:", paste(colnames(summary_df), collapse=", "), "\n")
