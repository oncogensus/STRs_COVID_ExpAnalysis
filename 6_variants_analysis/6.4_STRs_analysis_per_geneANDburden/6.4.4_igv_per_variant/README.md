# IGV.js por variante (navegador)

Workflow completo: gera BEDs + mapeamento BAM, depois sobe IGV.js para cada gene com outlier.

## Estrutura

```
6.4.4_igv_per_variant/
├── data/
│   ├── 0_generate_beds.R          # gera BEDs + TSV (R)
│   ├── 0_generate_beds.pbs        # PBS para rodar no cluster
│   ├── str_samples_bams.tsv       # output: mapeamento gene->BAM
│   ├── str_samples_with_variant.bed
│   └── str_samples_without_variant.bed
├── igv_variant.sh                 # IGV.js para 1 gene
├── run_all.sh                     # IGV.js para todos os genes
└── README.md
```

## Pipeline

```
suggestive_strs_outliers.tsv (6.4.1.2)
    ↓
data/0_generate_beds.R  →  data/*.bed + data/str_samples_bams.tsv
    ↓
igv_variant.sh GENE  →  IGV.js via HTTP
```

## 1. Gerar BEDs (no cluster)

```bash
cd 6.4.4_igv_per_variant/data
qsub 0_generate_beds.pbs
```

Ou localmente (se BAM dir acessivel):
```bash
Rscript 0_generate_beds.R
```

## 2. Rodar IGV.js — todos os genes

```bash
cd 6.4.4_igv_per_variant
bash run_all.sh
```

## 3. Rodar IGV.js — 1 gene

```bash
bash igv_variant.sh KCNQ5
bash igv_variant.sh KCNQ5 9000
```

## No PC (PowerShell)

```powershell
ssh -L 8201-82XX:localhost:8201-82XX Carlos_Chagas
```

Abra no navegador: `http://localhost:8201/tmp/igvjs_GENE/index.html`

## Pré-requisitos
- env `igv` no cluster (com `samtools`, `R`, `python`)
- BAM dir: `/storage/users/tulio/Projeto_Luy_COVID/results/recal/`
- Internet no PC para carregar igv.js do CDN

## Notas
- Os scripts são data-driven: leem `data/str_samples_bams.tsv` em runtime.
- BAMs sao extraidos por regiao (+/- 1000 bp) — nao servem BAM inteiro.
- Para remover arquivos temporarios: `rm -rf /tmp/igvjs_*`.
