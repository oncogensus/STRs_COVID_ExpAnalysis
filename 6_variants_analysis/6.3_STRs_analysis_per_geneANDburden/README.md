# 6.3 — STRs Analysis per Gene & Burden

Módulo principal de análise de STRs associados a COVID-19: cruzamento com DEGs de RNA-seq, visualização e inspeção via IGV.

## Estrutura

```
6.3_STRs_analysis_per_geneANDburden/
├── .gitignore
├── 6.3.1_pre_processing/          # (reservado para pré-processamento futuro)
├── 6.3.2_RNA_data_analysis/       # Cruzamento DEGs x STRs + análise estatística
│   ├── 6.3.2.2_RNA_matrix/
│   │   ├── 1_cross_intervention_STRs.py   # DEGs x STRs por intervenção
│   │   ├── 2_plot_intervention_summary.R  # Raincloud + tabela publicação
│   │   ├── 3_plot_raincloud_per_locus.R   # Raincloud por locus (outliers)
│   │   ├── 4_mann_whitney_allele.R        # Mann-Whitney U (case vs control)
│   │   ├── 5_no_overlap_analysis.R        # Análise de não-overlap
│   │   └── 8_outlier_detail_table.R       # Tabela detalhada de outliers
│   └── README.md
│
└── 6.3.4_igv_per_variant/          # Visualização IGV por variante
    ├── 1_generate_beds.R           # Gera BEDs + mapeamento BAM
    ├── 2_igv_variant.sh            # IGV.js para 1 STR
    ├── 3_run_all.sh                # Gera scripts IGV para todos os STRs
    ├── scripts/                    # Scripts IGV gerados por variante
    └── README.md
```

## Fluxo de Dados

```
samples/STRs_analysis_dataset.tsv (dataset unificado, etapa 6.1)
    ↓
1_cross_intervention_STRs.py  →  intervention_strs.tsv (todos os STRs em genes DEGs)
    ↓                            intervention_outliers.tsv (apenas STRs com outliers DBSCAN)
    ↓                            intervention_summary.tsv (resumo por gene/intervenção)
    ↓
2_plot_intervention_summary.R →  intervention_raincloud.png + publicação tables
    ↓
3_plot_raincloud_per_locus.R  →  Raincloud por locus (outliers DBSCAN apenas)
    ↓
4_mann_whitney_allele.R       →  Mann-Whitney U (case vs control)
    ↓
5_no_overlap_analysis.R       →  Loci sem overlap case/control
    ↓
8_outlier_detail_table.R      →  Tabela publication-ready detalhada
    ↓
1_generate_beds.R             →  BEDs + BAM mapping
    ↓
2_igv_variant.sh              →  IGV.js para inspeção visual
```

## Pré-requisitos

- `6.1_merge_datasets/` deve ser executado primeiro (gera `STRs_analysis_dataset.tsv`)
- `5_global_dbscan/` deve ter sido executado (outliers DBSCAN global)
- GSE subdiretórios com TSVs de DEGs (para 6.3.2)

## Ambiente

| Script | Ambiente (micromamba) |
|---|---|
| `1_cross_intervention_STRs.py` | `str` |
| `2_plot_intervention_summary.R` | `r_enrich_env` |
| `3_plot_raincloud_per_locus.R` | `r_enrich_env` |
| `4_mann_whitney_allele.R` | `r_enrich_env` |
| `5_no_overlap_analysis.R` | `r_enrich_env` |
| `8_outlier_detail_table.R` | `r_enrich_env` |
| `1_generate_beds.R` | `igv` |

## Execução

```bash
# 6.3.2: Cruzamento DEGs x STRs
cd 6.3.2_RNA_data_analysis/6.3.2.2_RNA_matrix
qsub 1_cross_intervention_STRs.pbs
# aguardar conclusão
qsub 2_plot_intervention_summary.pbs
qsub 3_submit_raincloud.pbs
qsub 4_mann_whitney_allele.pbs
qsub 5_no_overlap_analysis.pbs
qsub 8_outlier_detail_table.pbs

# 6.3.4: Visualização IGV
cd ../../6.3.4_igv_per_variant
qsub 1_generate_beds.pbs
bash 3_run_all.sh
```
