# UCSC Tracks Table

Builds the publication-ready UCSC track annotation table for the STR outlier loci associated with COVID-19 genes.

## Purpose

Combine the outlier detail table (from `6_outlier_detail_table.R`) with the structured query table (`main_table_structured.csv`) and enrich it with motif/purity (TRExplorer), RepeatMasker, ENCODE signals (DNase, ATAC, H3K4me3, H3K27ac, CTCF), JARVIS, cCRE and GeneHancer annotations.

## Structure

```
UCSC_tracks_table/
├── build_main_table.py       # Parser + table builder
├── main_table_structured.csv # Curated base table (input)
├── main_table_with_variants.csv   # Output: flat table
├── main_table_with_variants.html  # Output: gt publication table
├── main_table_with_variants.xlsx  # Output: Excel workbook
└── raw_data/                 # Per-gene source tables (input)
```

## Pipeline

```
main_table_structured.csv ─┐
raw_data/*.txt ─────────────┤→ build_main_table.py → main_table_with_variants.{csv,html,xlsx}
outlier_detail_table.html ──┘    (gene → region/chr/start/motif)
```

`build_main_table.py`:

1. Parses `outlier_detail_table.html` (default path; override with `--outlier-html`) to recover the genomic location of each gene.
2. Reads `main_table_structured.csv` (UTF-8 with cp1252/latin-1 fallback) and `raw_data/*.txt`.
3. Normalizes values: PT tissue names translated (via the `PT_TO_EN` dictionary), ENCODE `Max (Min)` → `min - max`, tissue-specific values merged.
4. Renders the final `great_tables` (gt) HTML table with grouped columns and source notes.

## Execution

```bash
# In the 6.3.2.2_RNA_matrix directory, after 6_outlier_detail_table.R:
python UCSC_tracks_table/build_main_table.py \
  --outlier-html results/outlier_detail_table.html
```

## Inputs

| File | Description |
|---|---|
| `outlier_detail_table.html` | Outlier detail table from `6_outlier_detail_table.R` |
| `main_table_structured.csv` | Curated per-locus query table |
| `raw_data/*.txt` | Per-gene source data (CAMK4, DRAIC, GNG7, MYOZ2, ROBO2, ST6GALNAC3) |

## Outputs

| File | Description |
|---|---|
| `main_table_with_variants.csv` | Flat enriched table |
| `main_table_with_variants.html` | Publication-ready `gt` table |
| `main_table_with_variants.xlsx` | Excel workbook with adjusted column widths |

## Notes

- `PT_TO_EN` maps Portuguese tissue names from the input tables to English (Lung / Brain / Lung and Brain).
- GeneHancer ID/score entries are curated per gene (currently `GNG7`).
- The table is labeled for the UCSC Genome Browser track set used in the publication.

## Environment

- Python with `beautifulsoup4`, `pandas`, `great_tables`, `openpyxl`