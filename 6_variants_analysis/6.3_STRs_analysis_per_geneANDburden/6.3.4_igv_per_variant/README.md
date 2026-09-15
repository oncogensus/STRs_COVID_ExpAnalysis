# IGV.js por variante (navegador)

Workflow completo: gera BEDs + mapeamento BAM, depois sobe IGV.js para cada STR com outlier.

## Estrutura

```
6.3.4_igv_per_variant/
├── 1_generate_beds.R          # gera BEDs + TSV (R)
├── 1_generate_beds.pbs        # PBS para rodar no cluster
├── str_samples_bams.tsv       # output: mapeamento STR->BAM
├── str_samples_with_variant.bed
├── str_samples_without_variant.bed
├── 2_igv_variant.sh             # IGV.js para 1 STR
├── 3_run_all.sh                 # IGV.js para todos os STRs
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

## 1. Gerar BEDs (no cluster)

```bash
cd 6.3.4_igv_per_variant
qsub 1_generate_beds.pbs
```

Ou localmente (se BAM dir acessivel):
```bash
Rscript 1_generate_beds.R
```

## 2. Rodar IGV.js — todos os STRs

```bash
cd 6.3.4_igv_per_variant
bash 3_run_all.sh
```

## 3. Rodar IGV.js — 1 STR

```bash
bash 2_igv_variant.sh chr1:76143392:GT:16
bash 2_igv_variant.sh chr1:76143392:GT:16 9000
```

## No PC (PowerShell)

```powershell
ssh -L 8201-82XX:localhost:8201-82XX Carlos_Chagas
```

Abra no navegador: `http://localhost:8201/tmp/igvjs_chr1_76143392_GT_16/index.html`

## Pré-requisitos
- env `igv` no cluster (com `samtools`, `R`, `python`)
- BAM dir: `/storage/users/tulio/Projeto_Luy_COVID/results/recal/`
- Internet no PC para carregar igv.js do CDN

## Notas
- Os scripts são data-driven: leem `str_samples_bams.tsv` em runtime.
- BAMs sao extraidos por regiao (+/- 1000 bp) — nao servem BAM inteiro.
- STRs_ID sao sanitizados para nomes de arquivo (`:` → `_`).
- Para remover arquivos temporarios: `rm -rf /tmp/igvjs_*`.
