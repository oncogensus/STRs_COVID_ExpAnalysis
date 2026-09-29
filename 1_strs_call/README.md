# 1 — STR Calling

Short Tandem Repeat (STR) genotyping of the cohort BAM files using [STRling](https://github.com/quinlan-lab/STRling-nf).

## Purpose

Extract STR-containing regions from the aligned BAMs, merge the per-sample bins, and call expanded STR loci.

## STRling Commands

```bash
strling extract -b sample.bam -f reference.fa -o output_dir
```

- `extract`: Extracts STR-containing regions from BAM files
- `-b sample.bam`: Input BAM file for analysis
- `-f reference.fa`: Reference genome in FASTA format
- `-o output_dir`: Output directory for results

## Required Files

- `hg38.fa`: Reference genome (FASTA)
- `hg38.fa.str`: Reference genome with STR metadata
- Sequenced sample files (BAM format)

## STRling Pipeline Stages

1. **Extract**: Collects STR regions from input BAM
2. **Merge**: Combines data from multiple samples (if joint analysis)
3. **Estimate**: Estimates STR expansions
4. **Call**: Calls expanded STR loci

## Structure

```
1_strs_call/
└── 1.1_strling_STRs_call.sh     # STRling extract, merge and call
```

## Configuration

| Variable | Default | Description |
|---|---|---|
| `PREFIX` | `*` | Glob filter for sample BAMs |
| `PREFIX_PATH` | `/storage/users/tulio/Projeto_Luy_COVID/results/recal` | Directory containing the aligned BAMs |
| `STRLING` | `/storage2/matheusbomfim/.../bin/strling` | Path to the STRling binary |
| `REF_FA` | `hg38.fa` | Reference genome (FASTA) |
| `GENOME_STR` | `hg38.fa.str` | Reference genome with STR metadata |
| `JOINT_DIR` | `cbgm_output` | Output directory for all STRling artifacts |

## Pipeline

1. **Extract** — collect STR regions from each `*.bam` into `<basefile>.str.bin`
2. **Merge** — combine all `.str.bin` files into a joint STRling data set
3. **Call** — call expanded STR loci per sample

The script aborts if no BAM matches the prefix or if any step fails.

## Execution

```bash
bash 1.1_strling_STRs_call.sh
```

## Outputs

- `cbgm_output/*.str.bin` — per-sample extracted STR bins
- `cbgm_output/strling-bounds.txt` — merged STR boundary file
- `cbgm_output/<sample>` — per-sample STR calling results (input for stage 2)

## Environment

- `str` (micromamba): STRling