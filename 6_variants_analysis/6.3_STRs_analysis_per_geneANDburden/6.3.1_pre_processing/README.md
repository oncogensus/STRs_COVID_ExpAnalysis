# COVID-19 Transcriptomic Analysis

Scripts used for preprocessing and differential expression analysis of
publicly available COVID-19 transcriptomic datasets from the Gene
Expression Omnibus (GEO).

## Datasets

### GSE157103 — Blood

Leukocyte RNA-seq data from COVID-19 patients.

Differential expression analyses included:

- COVID-19 ICU vs COVID-19 non-ICU
- HFD-45 = 0 vs HFD-45 > 0

Scripts for the dataset-specific preprocessing and differential expression
analysis are provided in `GSE157103/`.

### GSE183533 — Lung

Post-mortem lung transcriptomic data from severe COVID-19 cases and
uninfected controls.

Differential expression analysis:

- COVID-19 vs control

Scripts for the dataset-specific preprocessing and differential expression
analysis are provided in `GSE183533/`.

### GSE188847 — Brain

Post-mortem frontal cortex RNA-seq data from COVID-19 patients and
uninfected controls.

Differential expression analysis:

- COVID-19 vs control

Scripts for the dataset-specific preprocessing and differential expression
analysis are provided in `GSE188847/`.

## Repository structure

```text
6.3.1_pre_processing/
├── GSE157103/    # Blood
├── GSE183533/    # Lung
└── GSE188847/    # Brain
```

Each directory contains the scripts used for dataset-specific data preparation
and differential expression analysis.

Expression matrices, intermediate files, statistical results, figures, and
other generated outputs are not included in this repository.
