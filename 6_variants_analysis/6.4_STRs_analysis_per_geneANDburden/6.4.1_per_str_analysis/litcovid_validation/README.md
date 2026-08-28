# litcovid_validation  (Parte 1 da validação DBSCAN↔LitCovid)

Dados os STR loci marcados como **outliers pelo DBSCAN global** (coluna
`n_outliers > 0` no catálogo da coorte, produzido por `5_dbscan`), extraímos os
genes que os contêm e recuperamos a **literatura COVID do LitCovid** que menciona
cada gene. Responde a: "os genes que o sinal do DBSCAN aponta como relevantes têm
suporte na literatura COVID?"

## Etapas
1. `get_outlier_genes.py` — catálogo → `data/outlier_genes.txt` (gene, n_outlier_str_loci).
2. `download_litcovid.sh` — baixa `litcovid2pubtator.json.gz` (~2.3 GB, gene-anotado) → `data/`.
3. `gene_literature.py` — mapeia símbolo→Entrez (MyGene.info) e percorre o PubTator
   em streaming, coletando os artigos que mencionam cada gene-alvo →
   `results/gene_literature.tsv` + `results/gene_literature_summary.tsv`.

## Como rodar (cluster, env `igv`)
```bash
bash run_litcovid.sh               # ja inclui o download (1x, ~2.3GB) + get_outlier_genes + gene_literature
# ou via PBS: qsub litcovid.pbs  (submit_litcovid.sh)
```
Para submeter as DUAS etapas (Parte 1 -> Parte 2) encadeadas via PBS, use na pasta
pai: `bash submit_validation.sh`.
`CATALOG` default: `../../../samples/STRs_analysis_dataset.tsv` (sobrescreva se precisar).

## Saídas
- `results/gene_literature.tsv` — `gene, entrez, pmid, title, journal, year`.
- `results/gene_literature_summary.tsv` — `gene_name, n_outlier_str_loci, n_articles`.

## Dependências
- `python3` (micromamba `igv`); acesso à internet (MyGene.info + FTP NCBI).
- O catálogo completo de STRs da coorte (ver `CATALOG`).
