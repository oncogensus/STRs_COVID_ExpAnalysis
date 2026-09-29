# STR Analysis Workflow Development

![Python](https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)
![Status](https://img.shields.io/badge/Status-Development-yellow)

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Pipeline Architecture](#pipeline-architecture)
- [Module Documentation](#module-documentation)
- [Important Notes](#important-notes)
- [Stage 1: STR Calling](#stage-1-str-calling-strs_call)
- [Stage 2: Data Stratification](#stage-2-data-stratification-data_split)
- [Stage 3: Genomic Annotation](#stage-3-genomic-annotation-gtf_annot)
- [Stage 4: Ancestry Assignment](#stage-4-ancestry-assignment-ancestry)
- [Stage 5: Global DBSCAN Analysis](#stage-5-global-dbscan-analysis-5_global_dbscan)
- [Stage 6: Variant Analysis](#stage-6-variant-analysis-6_variants_analysis)
- [Workflow Execution](#workflow-execution)
- [Output Structure](#output-structure)
- [References](#references)

## Overview

This repository documents the **development workflow** for Short Tandem Repeat (STR) analysis throughout 2025-2026. The pipeline identifies STRs potentially impacting clinical outcomes in COVID-19 patients by comparing two groups without comorbidities:

- **Controls**: Survivors of severe COVID-19 (n=136, age <60)
- **Cases**: Fatal COVID-19 outcomes (n=32, age <60)

## Background

### Dataset Information

This analysis uses COVID-19 sequencing data from "Rare genetic variants and severe COVID-19 in previously healthy admixed Latin American adults" ([doi.org/10.1038/s41598-025-08416-1](https://doi.org/10.1038/s41598-025-08416-1)).

#### Patient Cohorts

- **Group 1** (Cases): Death without comorbidity, age < 60 years (n=32, mean age: 45.5 ± 19-59)
- **Group 2** (Controls): Survivors without comorbidity, severe COVID, age < 60 years (n=136, mean age: 43.1 ± 20-60)

### Tools and References

- [STRling](https://github.com/quinlan-lab/STRling-nf) — STR genotyping from short-read data ([documentation](https://strling.readthedocs.io/en/latest/index.html))
- [EthSEQ](https://github.com/mhguo1/AD_STR/tree/main) — Ancestry inference from sequencing data


---

## Pipeline Architecture

The analysis pipeline is organized as follows:

```text
STRs_COVID_ExpAnalysis/
├── 1_strs_call/
├── 2_data_split/
├── 3_gtf_annot/
├── 4_ancestry/
├── 5_global_dbscan/
└── 6_variants_analysis/
    ├── 6.1_merge_datasets/
    ├── 6.2_desc_data_viz/
    ├── 6.3_STRs_analysis_per_geneANDburden/
    │   ├── .gitignore
    │   ├── 6.3.1_pre_processing/
    │   ├── 6.3.2_RNA_data_analysis/
    │   │   └── 6.3.2.2_RNA_matrix/
    │   └── 6.3.4_igv_per_variant/
    └── 6.5_ancestry_analysis/
```

---

## Module Documentation

Each pipeline module has its own README with detailed inputs, outputs, and execution instructions:

| Module | Description | Documentation |
|---|---|---|
| 1 — STR Calling | STRling extract/merge/call | [`1_strs_call/README.md`](1_strs_call/README.md) |
| 2 — Data Stratification | Unify + quality filtering + case/control split | [`2_data_split/README.md`](2_data_split/README.md) |
| 3 — Genomic Annotation | GTF-based region + gene annotation | [`3_gtf_annot/README.md`](3_gtf_annot/README.md) |
| 4 — Ancestry Assignment | EthSEQ ancestry inference | [`4_ancestry/README.md`](4_ancestry/README.md) |
| 5 — Global DBSCAN | Normalization + DBSCAN outlier detection | [`5_global_dbscan/README.md`](5_global_dbscan/README.md) |
| 6.1 — Merge Datasets | Unified STR dataset | [`6.1_merge_datasets/README.md`](6_variants_analysis/6.1_merge_datasets/README.md) |
| 6.2 — Descriptive Analysis & Visualization | Coverage, DBSCAN validation, genome viz | [`6.2_desc_data_viz/README.md`](6_variants_analysis/6.2_desc_data_viz/README.md) |
| 6.3 — STRs per Gene & Burden | DEG x STR, raincloud plots, IGV | [`6.3_.../README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/README.md) |
| 6.5 — Ancestry Analysis | Ancestry x DBSCAN correlations + viz | [`6.5_ancestry_analysis/README.md`](6_variants_analysis/6.5_ancestry_analysis/README.md) |
| UCSC Tracks Table | Publication-ready table builder | [`.../UCSC_tracks_table/README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.3.2_RNA_data_analysis/6.3.2.2_RNA_matrix/UCSC_tracks_table/README.md) |

> The 6.3 module bundles sub-modules with their own READMEs: [`6.3.1_pre_processing/README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.3.1_pre_processing/README.md), [`6.3.2_RNA_data_analysis/README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.3.2_RNA_data_analysis/README.md), and [`6.3.4_igv_per_variant/README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/6.3.4_igv_per_variant/README.md).

---

## Stage 1: STR Calling (`strs_call`)

Identification of Short Tandem Repeat (STR) regions using STRling (extract, merge, call). See [`1_strs_call/README.md`](1_strs_call/README.md) for commands, required files, and pipeline stages.

---

## Stage 2: Data Stratification (`data_split`)

Separate identified STRs into control and case groups based on phenotype, with quality filtering (depth > 15, clips > 0, no homopolymers). See [`2_data_split/README.md`](2_data_split/README.md).

---

## Stage 3: Genomic Annotation (`gtf_annot`)

Annotation of genomic regions intersected by STRs (CDS, UTRs, promoters, introns, genes, intergenic), plus validation of the `others` category. See [`3_gtf_annot/README.md`](3_gtf_annot/README.md).

---

## Stage 4: Ancestry Assignment (`ancestry`)

Identification of global ancestry using EthSEQ (3D model). See [`4_ancestry/README.md`](4_ancestry/README.md).

---

## Stage 5: Global DBSCAN Analysis (`5_global_dbscan`)

Normalization of STR alleles by covariates (regression residuals) and per-locus outlier detection via DBSCAN. See [`5_global_dbscan/README.md`](5_global_dbscan/README.md).

---

## Stage 6: Variant Analysis (`6_variants_analysis`)

Integrated description, visualization, and filtering of identified variants.

- **6.1 Dataset Integration** — unified dataset. See [`6.1_merge_datasets/README.md`](6_variants_analysis/6.1_merge_datasets/README.md)
- **6.2 Descriptive Analysis & Visualization** — coverage, DBSCAN validation, genome viz. See [`6.2_desc_data_viz/README.md`](6_variants_analysis/6.2_desc_data_viz/README.md)
- **6.3 Per-STR Analysis** — DEG x STR, raincloud plots, IGV. See [`6.3_STRs_analysis_per_geneANDburden/README.md`](6_variants_analysis/6.3_STRs_analysis_per_geneANDburden/README.md)
- **6.5 Ancestry Analysis** — ancestry x DBSCAN correlations + viz. See [`6.5_ancestry_analysis/README.md`](6_variants_analysis/6.5_ancestry_analysis/README.md)

---

## Workflow Execution

### 1. Environment Setup

The pipeline runs in multiple conda/micromamba environments. `str_env.yaml`, `r_env.yaml` and `ethseq_env.yaml` are tracked in the repository (currently empty placeholders pending their final definitions); `r_viz`, `r_enrich_env`, `dbscan-r` and `igv` are pinned with the environments used on the cluster:

```bash
micromamba create -n str -f str_env.yaml
micromamba create -n r_env -f r_env.yaml
micromamba create -n ethseq_vcf_run -f ethseq_env.yaml
micromamba create -f r_viz.yaml
micromamba create -f r_enrich_env.yaml
micromamba create -f dbscan-r.yaml
micromamba create -f igv.yaml
```

See each module README for the environment required per script.

### 2. Prepare Reference Files

- Reference genome: `hg38.fa`
- Reference STR annotations: `hg38.fa.str`
- GTF annotations: `Homo_sapiens.GRCh38.98.gtf`
- Genome file:
  ```bash
  grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s' hg38.fa.fai | cut -f1,2 > genome.txt
  ```

### 3. Prepare Input Data

- BAM files for all samples
- Sample grouping file: `grupos.csv`
- Phenotype file: `samples_infos.csv`

### 4. Execution Order

1. STR Calling (`1_strs_call`)
2. Data Stratification (`2_data_split`)
3. Genomic Annotation (`3_gtf_annot`)
4. Ancestry Assignment (`4_ancestry`)
5. Global DBSCAN Analysis (`5_global_dbscan`)
6. Variant Analysis (`6_variants_analysis`)
   - 6.1 Dataset Integration (`6.1_merge_datasets/`)
   - 6.2 Descriptive Analysis & Visualization (`6.2_desc_data_viz/`)
   - 6.3 Per-STR Analysis (`6.3_STRs_analysis_per_geneANDburden/`)
   - 6.5 Ancestry Analysis (`6.5_ancestry_analysis/`)

## Output Structure

<details>
<summary>Click to expand full output tree</summary>

```text
STRs_COVID_ExpAnalysis/
├── 1_strs_call/
├── 2_data_split/
├── 3_gtf_annot/                    # Intermediate outputs from stages 1-3
│   └── samples/
│       ├── global_STRs_filtered.tsv
│       ├── summary_by_patient.tsv
│       ├── summary_report_final.tsv
│       ├── STRs_annotated_region.tsv
│       ├── global_annotation_statistics.tsv
│       └── others_regions_statistics.csv
├── 4_ancestry/                     # Ancestry outputs (stage 4)
│   └── EthSEQ_Results_3D/
│       └── Report.txt
├── 5_global_dbscan/                # DBSCAN results (stage 5)
│   ├── norm_test/
│   │   └── STRs_normalized_residuals.tsv
│   └── outliers_search/
│       └── results_dbscan/
│           └── outliers_per_str.tsv
└── 6_variants_analysis/            # Variant analysis (stage 6)
    ├── STRs_analysis_dataset.tsv   # Master integrated dataset
    ├── 6.1_merge_datasets/
    ├── 6.2_desc_data_viz/
    ├── 6.3_STRs_analysis_per_geneANDburden/
    │   ├── .gitignore
    │   ├── 6.3.1_pre_processing/
    │   ├── 6.3.2_RNA_data_analysis/
    │   │   └── 6.3.2.2_RNA_matrix/
    │   │       └── results/         # intervention_strs.tsv, intervention_outliers.tsv
    │   ├── 6.3.4_igv_per_variant/
    │   │   ├── str_samples_bams.tsv
    │   │   ├── str_samples_with_variant.bed
    │   │   ├── str_samples_without_variant.bed
    │   │   └── scripts/             # per-gene IGV scripts
    └── 6.5_ancestry_analysis/
        └── results/
            ├── categorical_data/
            ├── high_resolution/
            └── dataviz/
```

</details>

## References

1. STRling: Quinlan et al., *STRling: Rapid and accurate genotyping of short tandem repeats from short-read sequencing data*
2. EthSEQ: Liu et al., *EthSEQ: A tool for estimating ancestry from sequencing data*
3. COVID-19 HG: [https://www.covid19hg.org](https://www.covid19hg.org)
4. AD_STR DBSCAN Approach: [https://github.com/mhguo1/AD_STR](https://github.com/mhguo1/AD_STR/tree/main)

---

**Last Updated**: September 2026

**Status**: In development — planned 2026 publication
