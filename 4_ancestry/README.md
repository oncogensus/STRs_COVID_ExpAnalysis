# 4 — Ancestry Assignment

Infers global ancestry from the cohort genotypes using [EthSEQ](https://github.com/mhguo1/AD_STR/tree/main).

## Purpose

Estimate ancestry proportions (continuum) and likely population for each sample, used downstream for STR normalization (stage 5) and ancestry analyses (stage 6.5).

## Structure

```
4_ancestry/
├── 4.1_ethseq_vcf_run.r      # EthSEQ run (3D model)
└── EthSEQ_Results_3D/        # Output directory
    ├── Report.txt            # Per-sample ancestry assignment
    └── Report.PCAcoord       # PCA coordinates per sample
```

## Configuration

| Variable | Default | Description |
|---|---|---|
| `vcf_file` | — | Input gVCF/VCF to analyze |
| `model_available` | `Gencode.Exome` | EthSEQ model (see `getModelsList()`) |
| `model_assembly` | `hg38` | Reference genome |
| `model_pop` | `All` | Populations to evaluate |
| `out_dir` | `EthSEQ_Results_3D` | Output directory |
| `space` | `3D` | Ancestry space configuration |

## Execution

```bash
Rscript 4.1_ethseq_vcf_run.r
```

## Inputs

- Genotyped VCF of the cohort samples

## Outputs

- `EthSEQ_Results_3D/Report.txt` — ancestry assignment per sample (`pop`, `contribution`, `type`)
- `EthSEQ_Results_3D/Report.PCAcoord` — PCA coordinates used by the normalization step

## Environment

- `ethseq_vcf_run` (micromamba): `R`, `EthSEQ`