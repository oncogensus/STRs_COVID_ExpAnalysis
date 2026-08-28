# IGV.js por variante (navegador)

Página IGV.js para cada variante STR, sem necessidade de VNC.
Cada variante mostra:
- **anotação** apenas da amostra **com** a variante (`str_samples_with_variant.bed`, restrito a este locus)
- **reads (BAMs)** das **duas** amostras: variante + controle

## Pré-requisitos (no cluster)
- env `igv` ativo (ou `igv`/`python` no PATH)
- ter gerado `str_samples_bams.tsv` e `str_samples_with_variant.bed`:
  ```bash
  micromamba activate igv
  Rscript str_samples_to_bed.R
  ```
- `python` 3.7+ (suporte a Range, necessário para os BAMs)
- `samtools` (extrai a região do locus de cada BAM); instale se faltar:
  `micromamba install -n igv -c bioconda samtools`
- internet no PC para carregar o igv.js do CDN

## O que o script faz com os BAMs
Em vez de servir o BAM inteiro (que pode ter dezenas de GB), o script extrai
só a região do locus com `samtools view -b -h <bam> <chr>:inicio-fim` (flanco
de `+/-1000 bp`, ajustável pela variável `FLANK` no topo de cada script) e
reindexa. O IGV.js carrega então um BAM pequeno e rápido. Se a extração vier
vazia, ele usa o BAM completo como fallback.

## Uso — uma variante
No cluster:
```bash
bash igv_per_variant/ROBO2.sh
```
No seu PC (PowerShell):
```powershell
ssh -L 8201:localhost:8201 Carlos_Chagas
```
Abra no navegador: **http://localhost:8201/tmp/igvjs_ROBO2/index.html**

## Uso — todas de uma vez
No cluster:
```bash
bash igv_per_variant/run_all.sh
```
No PC:
```powershell
ssh -L 8201-8208:localhost:8201-8208 Carlos_Chagas
```
Abra as páginas em abas: `http://localhost:8201/tmp/igvjs_ROBO2/index.html` (ROBO2), `8202`→`/tmp/igvjs_ANK3/index.html`, … `8208`→`/tmp/igvjs_ST6GALNAC3/index.html`.

## Mapeamento variante → porta
| gene | porta |
|------|-------|
| ROBO2 | 8201 |
| ANK3 | 8202 |
| CDH12 | 8203 |
| NKAIN2 | 8204 |
| SEMA6D | 8205 |
| KCNH1 | 8206 |
| KCNQ5 | 8207 |
| ST6GALNAC3 | 8208 |

## Notas
- Os scripts são estáticos, mas data-driven: leem `str_samples_bams.tsv` em runtime, então continuam corretos se os BAMs mudarem.
- Ao fechar o script (Ctrl+C) o servidor HTTP daquela variante é encerrado.
- Para remover os arquivos temporários: `rm -rf /tmp/igvjs_*`.
