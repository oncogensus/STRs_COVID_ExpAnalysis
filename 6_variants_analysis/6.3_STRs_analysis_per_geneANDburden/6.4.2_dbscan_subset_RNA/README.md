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

O pipeline `6.4.3_burden_test` foi parametrizado para rodar sobre **ambos os
datasets** (GWAS-filtrado e RNA) e inclui um **script único** de comparação
(`compare_gwas_rna.R`, PBS `compare_gwas_rna.pbs`).

### Burden test — `6.4.3_burden_test`

`burden_gwas.R` aceita argumentos:

```
Rscript burden_gwas.R --strategy gwas_burden   # default, usa suggestive_gene_strs.tsv
Rscript burden_gwas.R --strategy rna_burden    # usa rna_gene_strs.tsv (RNA)
Rscript burden_gwas.R --strategy rna_burden --background <arquivo> --out-dir <dir>
```

- `gwas_burden` → `results_gwas_burden/`
- `rna_burden` → `results_rna/`

PBS: `burden_gwas.pbs` (GWAS), `burden_rna.pbs` (RNA) e `burden_both.pbs` (os dois em sequência + comparativo).

### Comparativo GWAS × RNA (`6.4.4_pathway_crossvalidation/compare_gwas_rna.R`)

Script único que substitui o antigo `compare_burden_hits.*` (removido) e os
scripts antigos de cross-validation da `6.4.4` (removidos). Gera em
`6.4.4_pathway_crossvalidation/results_gwas_rna_comparison/`:

- `strategy_outlier_sets.tsv` — flags por STR (outlier GWAS-sig em p<5e-8, outlier RNA, STRs dos genes de burden GWAS/RNA)
- `outlier_genes_union.tsv` — união de genes, nº de STRs por estratégia e overlap do maior alelo por gene
- `burden_hits_union_{uncorrected,corrected}.tsv` e `burden_hits_overlap_{uncorrected,corrected}.tsv` — hits SKAT em sobreposição/união (mesma lógica do comparativo antigo)
- `outlier_x_burden_by_gene.tsv` — matriz 2×2 outlier × burden-hit por gene
- `patient_str.tsv` — tabela longa STR × paciente (alelos `allele1_est`/`allele2_est`, `maior_alelo`, `group`, flags de fonte)
- `per_str_case_control.tsv` — descritivo por STR caso × controle (n, média/mediana/min/max/sd do maior alelo por grupo, nº de pacientes-outlier por grupo e `overlap_maior_alealo_grupos`) — **sem testes estatísticos**

Entradas: P1 GWAS (`covid_suggestive_genes_with_outlier_STRs.tsv`), outliers RNA
(`results/rna_outlier_genes.tsv`), SKAT das duas estratégias e o catálogo
`samples/STRs_analysis_dataset.tsv`.

PBS prontos: `compare_gwas_rna.pbs`; `burden_both.pbs` roda as duas estratégias e o comparativo em sequência.