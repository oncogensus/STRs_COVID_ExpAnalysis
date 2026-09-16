# 6.1 — Merge Datasets

Integra anotação genômica, resultados DBSCAN global, ancestralidade e dados demográficos em um dataset unificado.

## Objetivo

Consolidar as saídas dos estágios anteriores (anotação GTF, DBSCAN global, EthSEQ, fenótipo) em um único arquivo `STRs_analysis_dataset.tsv` que alimenta toda a etapa 6.

## Estrutura

```
6.1_merge_datasets/
└── 1_merge_datasets.r      # Script de integração
```

## Inputs

| Arquivo | Origem | Descrição |
|---|---|---|
| `STRs_annotated_region.tsv` | `3_gtf_annot/samples/` | Anotação genômica por STR |
| `outliers_per_str.tsv` | `5_global_dbscan/outliers_search/results_dbscan/` | Métricas DBSCAN global por STR |
| `Report.txt` | `4_ancestry/EthSEQ_Results_3D/` | Atribuição de ancestralidade |
| `samples_infos.csv` | `samples/` | Dados demográficos (grupo, idade, sexo) |

## Output

`samples/STRs_analysis_dataset.tsv` — dataset unificado com colunas:

| Categoria | Colunas |
|---|---|
| IDs & Clínico | `STRs_ID`, `group`, `age`, `sex` |
| Métricas STR | `allele1_est`, `allele2_est`, `depth` |
| Contexto Genômico | `repeat_unit`, `gene_id`, `gene_name`, `region`, `chrom`, `start`, `end`, `sample_id` |
| Genética Populacional | `pop`, `contribution`, `type` |
| Outliers DBSCAN Global | `n_outliers_dbscan_global`, `outlier_samples_dbscan_global`, `outlier_residuals_dbscan_global`, `n_clusters_dbscan_global`, `noise_ratio_dbscan_global` |

## Execução

```bash
cd 6.1_merge_datasets
Rscript 1_merge_datasets.r
```

## Ambiente

- `r_enrich_env` (micromamba): `data.table`, `dplyr`, `stringr`
