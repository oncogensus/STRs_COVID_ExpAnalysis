# IGV.js per variant (browser)

Full workflow: generates BEDs + BAM mapping, then launches IGV.js for each STR with outliers.

## Structure

```
6.3.4_igv_per_variant/
├── 1_generate_beds.R          # generates BEDs + TSV (R)
├── 1_generate_beds.pbs        # PBS to run on the cluster
├── str_samples_bams.tsv       # output: STR->BAM mapping
├── str_samples_with_variant.bed
├── str_samples_without_variant.bed
├── 2_igv_variant.sh             # IGV.js for 1 STR
├── 3_run_all.sh                 # IGV.js for all STRs
└── README.md
```

## Pipeline

```
intervention_outliers.tsv (6.3.2.2_RNA_matrix/results/)
    ↓
1_generate_beds.R  →  *.bed + str_samples_bams.tsv
    ↓
2_igv_variant.sh STRS_ID  →  IGV.js via HTTP
```

## 1. Generate BEDs (on the cluster)

```bash
cd 6.3.4_igv_per_variant
qsub 1_generate_beds.pbs
```

Or locally (if the BAM dir is accessible):
```bash
Rscript 1_generate_beds.R
```

## 2. Run IGV.js — all STRs

```bash
cd 6.3.4_igv_per_variant
bash 3_run_all.sh
```

## 3. Run IGV.js — 1 STR

```bash
bash 2_igv_variant.sh chr1:76143392:GT:16
bash 2_igv_variant.sh chr1:76143392:GT:16 9000
```

## On the PC (PowerShell)

```powershell
ssh -L 8201-82XX:localhost:8201-82XX Carlos_Chagas
```

Open in the browser: `http://localhost:8201/tmp/igvjs_chr1_76143392_GT_16/index.html`

## Prerequisites
- `igv` env (`samtools`, `R`, `python`), pinned as `igv.yaml` (repo root)
- BAM dir: `/storage/users/tulio/Projeto_Luy_COVID/results/recal/`
- Internet access on the PC to load igv.js from the CDN

## Notes
- The scripts are data-driven: they read `str_samples_bams.tsv` at runtime.
- BAMs are extracted per region (+/- 1000 bp) — the whole BAM is not served.
- STRs_ID are sanitized for file names (`:` → `_`).
- To remove temporary files: `rm -rf /tmp/igvjs_*`.