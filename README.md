# STR Analysis Workflow Development

## Overview

This repository documents the **ongoing development workflow** for Short Tandem Repeat (STR) analysis throughout 2025. It synthesizes the current stages of STR identification, outlier detection, and annotation designed to support exploratory analysis and variant filtering. The pipeline seeks to identify STRs potentially impacting clinical outcomes in COVID-19 patients by comparing two key groups without comorbidities:

- **Controls**: Survivors of severe COVID-19 (n=138, age <60)
- **Cases**: Fatal COVID-19 outcomes (n=33, age <60)

*Note: This workflow is actively under development for a planned 2026 publication evaluating STR associations with COVID-19 mortality.*

## Background

### Dataset Information

This analysis is based on COVID-19 sequencing data, the same dataset used for "Rare genetic variants and severe COVID-19 in previously healthy admixed Latin American adults" article (https://doi.org/10.1038/s41598-025-08416-1)

#### Patient Cohorts

- **Group 1**: Death without comorbidity, age < 60 years (n=33, mean age: 45.5 ± 19-59)
- **Group 2**: Survivors without comorbidity, severe COVID, age < 60 years (n=138, mean age: 43.1 ± 20-60)
- **Group 3**: Elderly patients (≥60 years) with comorbidity, mild COVID (n=36, mean age: 67.6 ± 60-82)

### Resources & References

- [STRling GitHub Repository](https://github.com/quinlan-lab/STRling-nf)
- [STRling Documentation](https://strling.readthedocs.io/en/latest/index.html)

---

## Pipeline Architecture

The analysis pipeline is organized as follows:

```
#/home/matheusbomfim/projects/strs_paper
├── strs_call/
├── data_split/
├── outliers_z-score/
├── gtf_annot/
├── ancestry/
├── desc_analysis/
└── results/
```

---

## Stage 1: STR Calling (`strs_call`)

### Purpose
Identification of Short Tandem Repeat regions using STRling.

### Command
```bash
strling extract -b sample.bam -f reference.fa -o output_dir
```

### Parameters
- `extract`: Extracts STR-containing regions from BAM files
- `-b sample.bam`: Input BAM file for analysis
- `-f reference.fa`: Reference genome in FASTA format
- `-o output_dir`: Output directory for results

### Required Files
- `hg38.fa`: Reference genome (FASTA)
- `hg38.fa.str`: Reference genome with STR metadata
- Sequenced sample files (BAM format)

### Environment
- STRling installed (micromamba environment: `str`)

### STRling Pipeline Stages
1. **Extract**: Collects STR regions from input BAM
2. **Merge**: Combines data from multiple samples (if joint analysis)
3. **Estimate**: Estimates STR expansions
4. **Call**: Calls expanded STR loci

**Script**: `strling_STRs_call.sh`

---

## Stage 2: Data Stratification (`data_split`)

### Purpose
Separate identified STRs into control and case groups based on phenotype.

### Required Files
- `grupos.csv`: Sample-to-group assignment file
- STRs identified from Stage 1 (output files)

### Output
Stratified data in `../samples/..` directories

**Script**: `group_split.sh`

---

## Stage 3: Statistical Outlier Detection (`outliers_z-score`)

### Purpose
Modified Z-score calculation comparing control versus case groups for STR outlier detection.

### Required Files
- `grupos.csv`: Sample group assignments
- `input_dir`: Sample locations (`../samples/gp_global`)
  - Genotypes files (`--genotypes`)
  - Unplaced variants

### Output
- `../samples/gp_global/outliers_gp_global/`: Z-score outlier results

**Script**: `strling-outliers_perGroup.py`

---

## Stage 4: Genomic Annotation (`gtf_annot`)

This stage unifies STR variants into consolidated files and annotates genomic regions using GTF reference annotations.

### 4.1: Global STR Consolidation

**Purpose**: Convert individual `--genotypes.txt` files into a unified consolidated file.

**Script**: `gp_global_merge.sh`

**Rationale**: Data stored in TSV format since subsequent analyses only require variant location, not VCF conversion.

**Required Files**
- `dir_vcfs`: Directory containing STRling-identified STRs (`*-genotype.txt`)
- `output_files`: Output directory for unified genotypes

**Output**
- `../samples/gp_global/merged/str_consolidated.tsv`: Unified STR genotypes

**Summary Statistics Generated**
- `variant_count`: Total number of STRs per sample
- `mean_depth`: Average sequencing depth
- `mean_allele_diff`: Mean allele size difference
- `homozygous_count`: Count of homozygous STRs
- `heterozygous_count`: Count of heterozygous STRs

### 4.2: Outlier STR Consolidation

**Purpose**: Merge individual Z-score outlier results into a unified file.

**Script**: `outliers_merge.sh`

**Rationale**: Data stored in TSV format for consistency with downstream analyses.

**Required Files**
- `dir_vcfs`: Directory containing Z-score outlier results (`outliers_gp_global*.STRs.tsv`)
- `output_files`: Output directory for unified outlier data

**Output**
- `../samples/gp_global/merged/outliers_consolidated.tsv`: Unified outlier STRs

**Summary Statistics Generated**
- `variant_count`: Count of outlier STRs per sample
- `mean_depth`: Average sequencing depth for outliers
- `mean_allele_diff`: Mean allele size difference for outliers
- `homozygous_count`: Count of homozygous outlier STRs
- `heterozygous_count`: Count of heterozygous outlier STRs

### 4.3: Global STR Genomic Annotation

**Purpose**: Annotate genomic locations of all identified STRs using GTF reference.

**Script**: `gtf_annot_global.py`

**Required Files**
- `TSV_PATH`: Location of STRs to be annotated
- `GTF_PATH`: GTF file with GRCh38.98 annotations
- `GENOME_FILE`: Chromosome size file (generated below)
- `OUTDIR`: Output directory for annotated results

**Genome File Generation**
```bash
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s' hg38.fa.fai | cut -f1,2 > genome.txt
```

**Environment**
- STRling installed (micromamba environment: `str`)

### 4.4: Z-Score Outlier Genomic Annotation

**Purpose**: Annotate genomic locations of Z-score detected STR outliers.

**Script**: `gtf_annot_z-score.py`

**Required Files**
- `TSV_PATH`: Location of outlier STRs to be annotated
- `GTF_PATH`: GTF file with GRCh38.98 annotations
- `GENOME_FILE`: Chromosome size file (generated as in 4.3)
- `OUTDIR`: Output directory for annotated results

**Environment**
- (micromamba environment: `str`)


---

## Stage 5: Ancestry Assignment (`ancestry`)

### Purpose
Determine global ancestry composition for each sample using EthSEQ.

**Script**: `ethseq_vcf_run.r`

### Required Files
- `vcf_file`: Genomic VCF file for ancestry inference

### Parameters
- `model_available`: Ancestry model selection
- `model_assembly`: Reference genome assembly (GRCh38)
- `model_pop`: Populations to evaluate
- `out_dir`: Output directory for ancestry results

**Environment**
- EthSEQ installed (micromamba environment: `ethseq`)


---

## Stage 6: Descriptive Analysis (`desc_analysis`)

### Purpose
Generate visualizations and descriptive statistics of STR distribution across genome and patient groups.

### 6.1: Data Extraction for Visualization

**Purpose**: Extract and prepare data for descriptive analysis.

**Script**: `desc_analysis_submit.r`

**Required Files**
- `samples_groups`: File mapping samples to groups
- `df_strs`: Annotated STRs file (`STRs_annotated_complete.tsv`)

### 6.2: Statistical Summaries and Visualizations

**Purpose**: Generate tables and plots describing STR distribution by position, repeat unit, and patient group.

**Script**: `desc_analysis.ipynb` (Jupyter Notebook)

**Required Input Files**

- `grupos.csv`: Sample-to-group assignment

**Input Data Files**
- `stats_group_repeat_alleles_meanallele_extension_quartiles_sd.csv`: Descriptive statistics (mean, SD, quartiles) of allele sizes and STR expansions stratified by group (case/control) and repeat type

- `strs_by_gene_combo.csv`: STRs collapsed to gene level with counts stratified by functional groups

- `strs_by_locus_combo.csv`: Locus-level data (one row per STR locus × sample) containing chromosome, genomic region, repeat unit, estimated alleles, and coverage depth

- `strs_by_region_combo.csv`: STRs aggregated by genomic region type (CDS, intron, promoter, etc.) with counts separated by group (case/control)

- `strs_by_repeatunit_combo.csv`: STRs aggregated by repeat unit motif showing total counts and group-stratified counts (n_case, n_control) for each motif

---

## Workflow Execution

### Prerequisites



1. Environment setup:
   
   This workflow uses three separate environments defined in YAML files (`str_env.yaml`, `r_env.yaml`, `ethseq_env.yaml`).

   Install environments from YAML files:

   **STRling environment** (for STR calling and annotation):
   ```
   micromamba create -f str_env.yaml
   micromamba activate str
   ```

   **R environment** (for statistical analysis and visualization):
   ```
   micromamba create -f r_env.yaml
   micromamba activate r_env
   ```

   **EthSEQ environment** (for ancestry inference):
   ```
   micromamba create -f ethseq_env.yaml
   micromamba activate ethseq
   ```

2. Prepare reference files
   - Reference genome (hg38.fa)
   - Reference STR annotations (hg38.fa.str)
   - GTF annotations (GRCh38.98)

3. Prepare input data
   - BAM files for all samples
   - Sample grouping file (grupos.csv)

### Recommended Execution Order

1. **STR Calling** (`strs_call`)
2. **Data Stratification** (`data_split`)
3. **Outlier Detection** (`outliers_z-score`)
4. **STR Consolidation** (`gtf_annot` sections 4.1-4.2)
5. **Genomic Annotation** (`gtf_annot` sections 4.3-4.4)
6. **Ancestry Analysis** (`ancestry`)
7. **Descriptive Analysis** (`desc_analysis`)

---

## Output Structure

All results are saved to the `../results` directory with the following organization:

```
results/
├── consolidated_strs/
│   ├── str_consolidated.tsv
│   ├── outliers_consolidated.tsv
│   ├── strs_annotated_complete.tsv
│   └── outliers_annotated_complete.tsv
├── ancestry/
│   └── [ancestry assignments and plots]
├── descriptive_analysis/
│   ├── stats_group_repeat_alleles_meanallele_extension_quartiles_sd.csv
│   ├── strs_by_gene_combo.csv
│   ├── strs_by_locus_combo.csv
│   ├── strs_by_region_combo.csv
│   ├── strs_by_repeatunit_combo.csv
│   └── [visualization plots]
```

## Required Reference Files (not included)
- hg38.fa: Download from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
- Homo_sapiens.GRCh38.98.gtf: Download from [Ensembl](http://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/)


---

## Key Decisions and Rationale

### TSV vs VCF Format
Data is maintained in TSV (or text) format rather than converted to VCF because downstream analyses only require variant genomic coordinates. VCF conversion would introduce unnecessary complexity without functional benefit.

### Z-Score Modification
The Z-score calculation is modified to account for group-specific distributions (control vs case) rather than applying a global threshold, enabling more sensitive detection of group-specific STR expansions.

### Consolidated Annotation Strategy
STRs are consolidated into unified files before genomic annotation to ensure consistent annotations across all samples and facilitate downstream comparative analyses.

---

## Citation and Attribution

If you use this pipeline in your research, please cite:

- STRling: [Quinlan et al.](https://github.com/quinlan-lab/STRling-nf)
- EthSEQ: Ancestry inference tool as referenced in Stage 5

---

## Contact & Support

For questions regarding this pipeline implementation, please refer to the original Notion documentation or contact the project maintainers.

---

**Last Updated**: December 2025  
**Status**: Workflow in devolpment for 2026 publication