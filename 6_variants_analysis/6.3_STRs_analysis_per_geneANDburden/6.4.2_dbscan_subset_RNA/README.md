# 6.4.2 — STRs em genes DEGs (RNA-seq) × DBSCAN

Cruzamento entre genes **diferencialmente expressos (DEGs)** de datasets de
RNA-seq públicos (**GSE157103**, **GSE188847**, **GSE183533**) e o catálogo de
STRs da coorte (STRling), com anotação de outliers DBSCAN (global e subset GWAS).

## Sub-etapas

### 6.4.2.2 — Matriz RNA × STRs (`6.4.2.2_RNA_matrix/`)

**`cross_DEGs_STRs.py`** cruza DEGs (Significant=Yes) de cada subpasta `GSE*` com
o catálogo de STRs e anota métricas DBSCAN.

**Entradas**:
- `--deg-dir` — raiz com subpastas `GSE*/` (cada GSE tem seus tsv/csv de DEGs)
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--gwas-outliers` — `6.4.1_dbscan_subset_GWAS/6.4.1.2_dbscan_subset_GWAS/results/suggestive_strs_outliers.tsv`
- `--out-dir` — diretório de saída (padrão `results/`)

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `all_STRs_in_DEGs.tsv` | todos os STRs anotados em genes DEGs |
| `outlier_STRs_in_DEGs.tsv` | apenas STRs com outliers DBSCAN global |
| `rna_gene_strs.tsv` | pares gene×STR dos DEGs (análogo GWAS `suggestive_gene_strs.tsv`), com coluna `datasets` (GSEs de origem) |
| `rna_outlier_genes.tsv` | STRs outlier por gene/DEG, **análogo GWAS `covid_suggestive_genes_with_outlier_STRs.tsv`**, com `dataset` e `gse` de origem |
| `rna_outlier_genes_by_study.tsv` | contagem de STRs outlier por (gene, GSE) |
| `rna_summary_by_study.tsv` | resumo pós-estudo por (GSE, gene): STRs identificadas (+outliers) e se há sobreposição do alelo maior entre grupos (`sim`/`nao`/`sem_dados`) |

Cada registro das saídas RNA traz a origem do estudo: coluna `dataset`
(`GSE157103/<arquivo>`), `gse` (`GSE157103`) ou `datasets` (GSEs separados por `;`).

**Submissão**: `cross_DEGs_STRs.pbs`.

---

## Avaliações entre GWAS × RNA

Os pipelines `6.4.3_burden_test` e `6.4.4_pathway_crossvalidation` foram
parametrizados para rodar sobre **ambos os datasets** (GWAS-filtrado e RNA).

### Burden test — `6.4.3_burden_test`

`burden_gwas.R` agora aceita argumentos:

```
Rscript burden_gwas.R --strategy gwas_burden   # default, usa suggestive_gene_strs.tsv
Rscript burden_gwas.R --strategy rna_burden    # usa rna_gene_strs.tsv (RNA)
Rscript burden_gwas.R --strategy rna_burden --background <arquivo> --out-dir <dir>
```

- `gwas_burden` → `results_gwas_burden/`
- `rna_burden` → `results_rna/`

PBS: `burden_gwas.pbs` (GWAS) e `burden_rna.pbs` (RNA).

### Pathway / gene cross-validation — `6.4.4_pathway_crossvalidation`

Os scripts `.R` aceitam `--p1-file`, `--p2-file`, `--p1-label`, `--p2-label`,
`--out` (e `--gene-out`/`--fig` no plot). Comparação GWAS × RNA (pré-configurada):

```
Rscript 1_pathway_crossvalidation.R --p2-file .../results/rna_outlier_genes.tsv --p1-label GWAS --p2-label RNA --out pathway_convergence_gwas_rna
Rscript 2_gene_crossvalidation.R    --p2-file .../results/rna_outlier_genes.tsv --p1-label GWAS --p2-label RNA --out gene_convergence_gwas_rna
Rscript 3_plot_pathways.R           --p2-file .../results/rna_outlier_genes.tsv --p1-label GWAS --p2-label RNA --out pathway_convergence_gwas_rna --gene-out gene_convergence_gwas_rna --fig gwas_rna_
```

**Pré-requisito**: `cross_DEGs_STRs.py` já rodado (gera `results/rna_outlier_genes.tsv`).

PBS prontos: `1_pathway_crossvalidation_rna.pbs`, `2_gene_crossvalidation_rna.pbs`,
`3_plot_pathways_rna.pbs`, e o orquestrador `submit_pathway_crossvalidation_rna.sh`.
(Os `.pbs` originais, sem sufixo `_rna`, preservam o comportamento GWAS × agnóstico.)