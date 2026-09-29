# 3 — Genomic Annotation

Annotates each STR with its genomic region using a GTF-based hierarchy (CDS, UTRs, promoters, introns, genes, intergenic).

## Purpose

Classify every filtered STR into a genomic region and attach gene metadata (`gene_id`, `gene_name`, `gene_biotype`), with a dedicated debug script auditing the "others" category.

## Structure

```
3_gtf_annot/
├── 3.1_gtf_annot_global.py      # Main annotation pipeline
├── 3.2_gtf_annot_global_debug.py # Audit of "others" + integrity checks
└── genome.txt                   # Chromosome sizes (chr1-22, X, Y)
```

## 3.1 — Global STR Annotation

Hierarchical annotation by region priority:

| Priority | Region | Description |
|---|---|---|
| 1 | `CDS` | Coding sequence |
| 2 | `five_prime_utr` | 5' UTR |
| 3 | `three_prime_utr` | 3' UTR |
| 4 | `non_coding_exons` | Non-translated exons |
| 5 | `promoter` | Upstream of TSS (default 3 kb) |
| 6 | `intron` | Intronic regions |
| 7 | `others` | Non-coding genes |
| 8 | `intergenic` | Regions between genes |

### Execution

```bash
python 3.1_gtf_annot_global.py
```

### Inputs

| File | Source |
|---|---|
| `global_STRs_filtered.tsv` | Stage 2 output |
| `Homo_sapiens.GRCh38.98.gtf` | Ensembl GTF annotation |
| `genome.txt` | Chromosome sizes (see below) |

### Outputs

- `STRs_annotated_region.tsv` — STRs with genomic annotation
- `code_regions_statistics.tsv` — summary and distribution statistics

### Generating `genome.txt`

```bash
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s' hg38.fa.fai | cut -f1,2 > genome.txt
```

## 3.2 — Annotation Validation (debug)

Audits variants classified as `others`:

- Region distribution of all variants
- Biotype distribution within `others`
- Integrity check (`others`/`intergenic` must have `gene_id = .`)
- Sample of `others` variants
- Protein-coding overlap audit against the GTF (via `BedTool` intersection)

### Execution

```bash
python 3.2_gtf_annot_global_debug.py
```

### Outputs

- `others_regions_statistics.csv` — detailed table of `others` variants
- `temp_others.bed` / `temp_coding.bed` — temporary BED files (auto-removed)

## Environment

- `str` (micromamba): `polars`, `pybedtools`