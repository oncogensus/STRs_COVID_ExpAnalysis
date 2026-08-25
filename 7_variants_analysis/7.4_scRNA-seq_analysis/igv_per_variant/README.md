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
- internet no PC para carregar o igv.js do CDN

## Uso — uma variante
No cluster:
```bash
bash igv_per_variant/ROBO2.sh
```
No seu PC (PowerShell):
```powershell
ssh -L 8201:localhost:8201 Carlos_Chagas
```
Abra no navegador: **http://localhost:8201**

## Uso — todas de uma vez
No cluster:
```bash
bash igv_per_variant/run_all.sh
```
No PC:
```powershell
ssh -L 8201-8208:localhost:8201-8208 Carlos_Chagas
```
Abra as páginas em abas: `http://localhost:8201` (ROBO2), `8202` (ANK3), … `8208` (ST6GALNAC3).

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
- Para remover: `rm -rf igv_per_variant/data`.
