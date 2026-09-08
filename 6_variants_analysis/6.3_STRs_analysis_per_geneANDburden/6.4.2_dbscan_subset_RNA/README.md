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

### 6.4.2.3 — Visualização RNA × STRs (`6.4.2.3_rna_visualization/`)

**`plot_rna_summary.R`** gera visualizações para publicação a partir dos
resultados de `cross_DEGs_STRs.py` e do catálogo de STRs.

**Entradas**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv` (com coluna `group`: case/control)
- `--rna-gene-strs` — `rna_gene_strs.tsv` (pares gene×STR com `datasets`)
- `--rna-outliers` — `rna_outlier_genes.tsv` (outliers DBSCAN global)
- `--summary` — `rna_summary_by_study.tsv` (resumo por estudo×gene)
- `--out-dir` — diretório de saída (padrão `results/`)

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `rna_ridgeline_by_study.png` | Ridgeline plot: allele2_est × densidade, facet por GSE, cor = case/control, triângulos pretos = outliers DBSCAN |
| `rna_publication_table.tsv` | Tabela por gene: n_STRs, n_outliers, n_overlap, proporções |
| `rna_publication_table.html` | Tabela `gt` formatada para publicação |

**Submissão**: `plot_rna_summary.pbs`.

### 6.4.2.3b — Teste estatístico Mann-Whitney (`mann_whitney_allele.R`)

Teste de Mann-Whitney (Wilcoxon rank-sum) para comparar `allele2_est` entre
grupos case e control, para cada STR localizado em genes DEGs.

**Entradas**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv` (com `group`: case/control)
- `--rna-gene-strs` — `rna_gene_strs.tsv`
- `--out-dir` — diretório de saída (padrão `results/`)

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `mann_whitney_mean_allele_results.tsv` | Tabela mean_allele: U, p, p_adjusted (BH), effect_size_r, medias/grupos |
| `mann_whitney_allele2_results.tsv` | Tabela allele2: mesmas colunas |
| `mann_whitney_manhattan_mean_allele.png` | Manhattan plot mean_allele: -log10(p) por STR |
| `mann_whitney_manhattan_allele2.png` | Manhattan plot allele2: -log10(p) por STR |
| `mann_whitney_boxplot_mean_allele.png` | Boxplot mean_allele case vs control (significativos BH<0.05) |
| `mann_whitney_boxplot_allele2.png` | Boxplot allele2 case vs control (significativos BH<0.05) |
| `mann_whitney_concordance.tsv` | Comparativo entre metricas (concordancia/discordancia) |
| `mann_whitney_concordance_plot.png` | Scatter effect_size mean_allele vs allele2, colorido por concordancia |

**Submissão**: `mann_whitney_allele.pbs`.

### 6.4.2.3c — Radar: localização genômica (`radar_genomic_location.R`)

Radar plots mostrando a distribuição de outliers e variantes sem sobreposição
por **região genômica** (promoter, intron, UTR, intergenic, etc.), com facet
por GSE study + radar combinado.

**Entradas**:
- `--rna-outliers` — `rna_outlier_genes.tsv`
- `--summary` — `rna_summary_by_study.tsv`
- `--out-dir` — diretório de saída (padrão `results/`)

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `radar_outliers_by_study.png` | Radar outliers DBSCAN: facet por GSE + combinado (vermelho) |
| `radar_no_overlap_by_study.png` | Radar sem sobreposição: facet por GSE + combinado (amarelo) |
| `radar_genomic_summary.tsv` | Tabela: região × GSE × categoria (outliers/no_overlap) |

**Submissão**: `radar_genomic_location.pbs`.

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