# IGV.js por variante (navegador)

Script generico IGV.js para qualquer gene com outlier STR.
Le `str_samples_bams.tsv` em runtime — funciona para qualquer gene sem hardcoded.

## Pipeline

```
6.4.1.2 (DBSCAN subset) → suggestive_strs_outliers.tsv
    ↓
6.4.5_str_samples_to_bed.R → str_samples_bams.tsv + str_samples_with_variant.bed
    ↓
6.4.4_igv_per_variant/ → IGV.js para cada gene
```

## Gerar os dados (no cluster)

```bash
cd /storage2/matheusbomfim/projects/git_repos/STRs_COVID_Analysis
qsub 6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.5_str_samples_to_bed.pbs
```

Ou localmente (se BAM dir acessivel):
```bash
Rscript 6.4.5_str_samples_to_bed.R
```

## Rodar IGV.js — todos os genes

No cluster:
```bash
cd 6_variants_analysis/6.4_STRs_analysis_per_geneANDburden/6.4.4_igv_per_variant
bash run_all.sh
```

No PC (PowerShell):
```powershell
ssh -L 8201-82XX:localhost:8201-82XX Carlos_Chagas
```

Abra no navegador: `http://localhost:8201/tmp/igvjs_GENE/index.html`

## Rodar IGV.js — um gene

```bash
bash igv_variant.sh KCNQ5
bash igv_variant.sh KCNQ5 9000   # porta especifica
```

## O que o script faz com os BAMs
Em vez de servir o BAM inteiro, extrai so a regiao do locus com
`samtools view -b -h <bam> <chr>:inicio-fim` (flanco +/-1000 bp) e
reindexa. O IGV.js carrega um BAM pequeno e rapido.

## Notas
- Os scripts sao data-driven: leem `str_samples_bams.tsv` em runtime.
- Ao fechar o script (Ctrl+C) o servidor HTTP daquela variante e encerrado.
- Para remover arquivos temporarios: `rm -rf /tmp/igvjs_*`.
- BAM dir: `/storage/users/tulio/Projeto_Luy_COVID/results/recal/` (ajustavel em 6.4.5).
