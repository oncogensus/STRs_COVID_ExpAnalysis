# covid19hg_evaluation

Cruzamento entre os **genes significativos do COVID-19 Host Genetics Initiative (r7)**
e o **catálogo completo de STRs da coorte** (STRling).

## Objetivo
1. Extrair, de forma *data-driven*, todos os genes com associação
   genome-wide significativa no COVID-19 HG r7 (fenótipos A2 = infecção/susceptibilidade,
   B2 = hospitalização, C2 = caso crítico).
2. Verificar em quais desses genes a coorte tem STR(s).

## Dados
- **COVID-19 HG r7** ("full" `leave_23andme`, completos, sem 23andMe):
  `COVID19_HGI_{A2,B2,C2}_ALL_leave_23andme_20220403.tsv.gz`
  (base: `https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats/`).
  `#CHR` é numérico (`3`, `10`, `X=23`).
- **Genes hg38**: GENCODE v39 primary_assembly GTF.
- **Catálogo de STR da coorte**: `STRs_analysis_dataset.tsv`
  (saída STRling; colunas `chrom/left/right/repeatunit/sample`).
  Ajuste `CATALOG` em `run_all.sh`.

## Método
- SNPs com `p < 5e-8` são mapeados ao gene se caem no corpo do gene ou até
  **50 kb** flanqueando (configurável: `--window`).
- O catálogo de STR é deduplicado por locus (`chrom-left-right-repeatunit`) e
  cada STR é cruzado por coordenada com os intervalos de gene (gene ± 50 kb).
- Cruzamento coord-a-coord (sem depender de nome de gene).

## Scripts
| script | função |
|---|---|
| `download_data.sh` | baixa summary stats + GTF para `./data` |
| `build_gene_bed.py` | GTF → `genes.hg38.bed` (corpos gênicos) |
| `extract_covid_genes.py` | summary → `covid_genes.tsv` |
| `overlap_str_cohort.py` | `covid_genes.tsv` × catálogo → `covid_genes_with_cohort_STRs.tsv` |
| `run_all.sh` | orquestra tudo |

## Como rodar (no cluster, env `igv`)
```bash
bash download_data.sh            # ou pule se já tiver ./data
python3 build_gene_bed.py --gtf data/gencode.v39.primary_assembly.annotation.gtf.gz --out genes.hg38.bed
python3 extract_covid_genes.py --sumstats data/*.leave_23andme_20220403.tsv.gz \
        --gene-bed genes.hg38.bed --out covid_genes.tsv
python3 overlap_str_cohort.py --catalog "$CATALOG" --covid-genes covid_genes.tsv \
        --out covid_genes_with_cohort_STRs.tsv
# ou simplesmente: bash run_all.sh   (ajuste CATALOG antes)
```

## Verificação de sanidade
`covid_genes.tsv` deve conter genes conhecidos do COVID-19 (ex.: `ACE2`, `OAS1`,
`TMPRSS2`, `ABO`, `MUC5B`, `DPP9`, `TYK2`). Se faltarem, revise parse de colunas/p.

## Caveats
- Arquivos `leave_23andme` excluem a 23andMe (menos poder, mas são os
  summary stats completos disponíveis).
- Mapeamento SNP→gene por proximidade (±50 kb) pode atribuir a múltiplos genes
  em regiões densas; janela configurável.
- `STRs_analysis_dataset.tsv` é pós-QC do STRling (depth≥15, sem homopolímeros).
- Todos os builds são hg38; cromossomos normalizados para `chrN` internamente.

## PBS (job scheduler do cluster)
Cada etapa tem um wrapper `.pbs` (fila `workq`, env `igv` via micromamba).
Submissão em cadeia (com dependência `afterok`):

```bash
bash submit_all.sh
# ou, manualmente:
qsub download_data.pbs
qsub build_gene_bed.pbs
qsub extract_covid_genes.pbs
qsub overlap_str_cohort.pbs
```

`submit_all.sh` roda no nó de login e encadeia as 4 etapas. O download
(`download_data.pbs`) tem `walltime=12h`; ajuste se necessário. O caminho do
catálogo de STR em `overlap_str_cohort.pbs`/`run_all.sh` deve bater com o seu.
