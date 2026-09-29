# 2 — Data Stratification

Unifies the per-sample STR calling outputs and stratifies the variants into case/control groups with quality filters.

## Purpose

Merge all `-genotype.txt` files, attach patient metadata (group, sample), apply quality filtering, and produce summary tables for the downstream stages.

## Structure

```
2_data_split/
├── 2.1_strling_unify_workflow.py  # Unification + quality filtering
├── grupos.csv                     # Sample-to-group assignment (input)
├── str_workflow.out               # Execution log (stdout)
└── str_workflow.err               # Execution log (stderr)
```

## Quality Filters

- `--min_depth` — minimum sequencing depth (default: **15**)
- `--min_clips_sum` — minimum sum of left + right clips (default: **1**, i.e. > 0)
- `--exclude_homopolymers` — removes loci with a repeat unit of length 1

## Execution

```bash
python 2.1_strling_unify_workflow.py \
  --input_dir <STRling -genotype.txt dir> \
  --groups grupos.csv \
  --out_dir ../samples
```

## Inputs

| File | Description |
|---|---|
| `<sample>-genotype.txt` | STRling calling output for each sample |
| `grupos.csv` | Mapping of `sample` to `group` (case/control) |

## Outputs

| File | Description |
|---|---|
| `global_STRs_filtered.tsv` | Quality-filtered STRs (input for stage 3) |
| `summary_by_patient.tsv` | Per patient: variant count, mean depth, zygosity |
| `summary_report_final.tsv` | Per group + global summary metrics |

## Environment

- `str` (micromamba): `pandas`