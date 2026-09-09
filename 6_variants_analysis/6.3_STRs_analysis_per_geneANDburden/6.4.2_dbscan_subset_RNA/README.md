# 6.4.2 — STRs em genes DEGs (RNA-seq) × DBSCAN

Cruzamento entre genes **diferencialmente expressos (DEGs)** de datasets de
RNA-seq públicos (**GSE157103**, **GSE188847**, **GSE183533**) e o catálogo de
STRs da coorte (STRling), com anotação de outliers DBSCAN (global e subset GWAS).

Análises organizadas em duas categorias:
- **Per study (GSE)** — análise por estudo de origem
- **Per intervention** — análise por comparação/intervenção (ex: COVID_ICU_vs_NonICU)

## Estrutura

```
6.4.2_dbscan_subset_RNA/
├── 6.4.2.2_RNA_matrix/
│   ├── per_study/
│   │   ├── cross_DEGs_STRs.py          Cruzação DEGs × STRs (por GSE)
│   │   └── cross_DEGs_STRs.pbs
│   │
│   └── per_intervention/
│       ├── cross_intervention_STRs.py   Cruzação DEGs × STRs (por intervenção)
│       └── cross_intervention_STRs.pbs
│
└── 6.4.2.3_rna_visualization/
    ├── per_study/
    │   ├── plot_rna_summary.R           Raincloud + tabela publication (por GSE)
    │   ├── plot_rna_summary.pbs
    │   ├── mann_whitney_allele.R        Mann-Whitney U test
    │   ├── mann_whitney_allele.pbs
    │   ├── radar_genomic_location.R     Radar: localização genômica
    │   └── radar_genomic_location.pbs
    │
    └── per_intervention/
        ├── plot_intervention_summary.R   Raincloud + tabela publication (por intervenção)
        └── plot_intervention_summary.pbs
```

---

## 6.4.2.2 — Matriz RNA × STRs

### Per study: `cross_DEGs_STRs.py`

Cruza DEGs (Significant=Yes) de cada subpasta `GSE*` com o catálogo de STRs.

**Entradas**:
- `--deg-dir` — raiz com subpastas `GSE*/`
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--gwas-outliers` — `6.4.1.2_.../results/suggestive_strs_outliers.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `all_STRs_in_DEGs.tsv` | todos os STRs anotados em genes DEGs |
| `outlier_STRs_in_DEGs.tsv` | apenas STRs com outliers DBSCAN global |
| `rna_gene_strs.tsv` | pares gene×STR com `datasets` (GSEs de origem) |
| `rna_outlier_genes.tsv` | STRs outlier por gene/DEG com `dataset` e `gse` |
| `rna_outlier_genes_by_study.tsv` | contagem de STRs outlier por (gene, GSE) |
| `rna_summary_by_study.tsv` | resumo por (GSE, gene): STRs, outliers, overlap |

**Submissão**: `qsub per_study/cross_DEGs_STRs.pbs`

### Per intervention: `cross_intervention_STRs.py`

Cruza DEGs de cada intervenção (arquivo TSV por comparação) com o catálogo de STRs.
Nome da intervenção extraído do nome do arquivo (prefixo DEG(s)_ e sufixos removidos).

**Entradas**:
- `--deg-dir` — raiz com subpastas `GSE*/` (cada GSE tem TSVs de intervenções)
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `intervention_strs.tsv` | todos os STRs anotados por intervenção |
| `intervention_outliers.tsv` | STRs outliers por intervenção |
| `intervention_summary.tsv` | resumo por (intervenção, gene): STRs, outliers, overlap |

**Submissão**: `qsub per_intervention/cross_intervention_STRs.pbs`

---

## 6.4.2.3 — Visualização RNA × STRs

### Per study: `plot_rna_summary.R`

Raincloud plot (apenas outliers DBSCAN) + tabela publication ready por GSE.

**Entradas**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--rna-gene-strs` — `rna_gene_strs.tsv`
- `--rna-outliers` — `rna_outlier_genes.tsv`
- `--summary` — `rna_summary_by_study.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `rna_raincloud_by_study.png` | Raincloud: allele2_est por GSE, cor = case/control |
| `rna_publication_table.tsv` | Tabela por GSE (1 row/GSE) |
| `rna_publication_table.html` | Tabela `gt` publication ready (Arial, spanners) |

**Submissão**: `qsub per_study/plot_rna_summary.pbs`

### Per intervention: `plot_intervention_summary.R`

Raincloud plot + tabela publication ready por intervenção.

**Entradas**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--intervention-strs` — `intervention_strs.tsv`
- `--intervention-outliers` — `intervention_outliers.tsv`
- `--intervention-summary` — `intervention_summary.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `intervention_raincloud.png` | Raincloud: allele2_est por intervenção, cor = case/control |
| `intervention_publication_table.tsv` | Tabela por intervenção (1 row/intervenção) |
| `intervention_publication_table.html` | Tabela `gt` publication ready (Arial, spanners) |

**Submissão**: `qsub per_intervention/plot_intervention_summary.pbs`

### 6.4.2.3b — Teste estatístico Mann-Whitney (`per_study/mann_whitney_allele.R`)

Teste de Mann-Whitney (Wilcoxon rank-sum) para comparar `allele2_est` entre
grupos case e control, para cada STR localizado em genes DEGs.

**Entradas**:
- `--str-catalog` — `samples/STRs_analysis_dataset.tsv`
- `--rna-gene-strs` — `rna_gene_strs.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `mann_whitney_mean_allele_results.tsv` | Tabela mean_allele: U, p, p_adjusted (BH), effect_size_r |
| `mann_whitney_allele2_results.tsv` | Tabela allele2: mesmas colunas |
| `mann_whitney_manhattan_*.png` | Manhattan plots |
| `mann_whitney_boxplot_*.png` | Boxplots (significativos BH<0.05) |
| `mann_whitney_concordance.tsv` | Comparativo entre métricas |
| `mann_whitney_concordance_plot.png` | Scatter effect_size |

**Submissão**: `qsub per_study/mann_whitney_allele.pbs`

### 6.4.2.3c — Radar: localização genômica (`per_study/radar_genomic_location.R`)

Radar plots: distribuição de outliers e variantes sem sobreposição
por **região genômica**, facet por GSE + radar combinado.

**Entradas**:
- `--rna-outliers` — `rna_outlier_genes.tsv`
- `--summary` — `rna_summary_by_study.tsv`
- `--out-dir` — diretório de saída

**Saídas** em `--out-dir`:
| Arquivo | Conteúdo |
|---|---|
| `radar_outliers_by_study.png` | Radar outliers DBSCAN (vermelho) |
| `radar_no_overlap_by_study.png` | Radar sem sobreposição (amarelo) |
| `radar_genomic_summary.tsv` | Tabela: região × GSE × categoria |

**Submissão**: `qsub per_study/radar_genomic_location.pbs`

---

## Ordem de execução

### Pipeline per study (GSE):
```bash
git pull
qsub per_study/cross_DEGs_STRs.pbs
qsub per_study/plot_rna_summary.pbs
qsub per_study/mann_whitney_allele.pbs
qsub per_study/radar_genomic_location.pbs
```

### Pipeline per intervention:
```bash
git pull
qsub per_intervention/cross_intervention_STRs.pbs
qsub per_intervention/plot_intervention_summary.pbs
```
