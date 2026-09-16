# burden_test — Genetic Association (burden / SKAT) — Etapa 7.6 do pipeline

Teste de associação entre **status de outlier em STR** e **mortalidade por
COVID-19** (caso/controle binário: óbito vs. sobrevivente, ambos < 60 anos, sem
comorbidades), ajustando por `age`, `sex` e 3 componentes principais de
ancestralidade do EthSEQ (`EV1`–`EV3`).

## Duas estratégias (definição de "outlier")

| script | estratégia | fonte de outliers | `remove_sample_outliers` |
|---|---|---|---|
| `burden_association.R` | DBSCAN | `7.4.1_per_str_analysis/outliers_search/results_dbscan/outliers_per_str.tsv` | `TRUE` |
| `burden_association_gwas.R` | GWAS-based (COVID-19 HG r7) | `covid19hg_evaluation/dbscan_subset/results/suggestive_strs_outliers.tsv` | `FALSE` |

- **DBSCAN**: outliers no espaço de resíduos normalizados completo; amostras
  marcadas como outliers globais (>2 DP) são removidas.
- **GWAS-based**: apenas STRs em genes sugestivos do COVID-19 HG r7 (p < 1e-5)
  com outlier na coorte. Cada STR tem 1 outlier → remoção global desligada
  (`FALSE`) para não apagar o sinal.

## Métodos
- **Burden por gene**: nº de STRs outlier por amostra dentro do gene como
  preditor (logística; covariáveis idem). `min_strs_per_gene = 2` (exclusão).
- **SKAT** por gene (`SKAT::SKAT`, `linear.weighted`, `out_type="D"`),
  ajuste small-sample automático (n=168 < 2000).
- **Burden por STR** (`run_str_burden_test`): regressão logística univariada
  por STR → `str_burden.tsv`.
- Toda saída por gene traz coluna `str_ids` (`;`-separada de `STRs_ID`) para
  rastrear hits; STRs intergênicas têm `gene_name="."`.

> Construção da matriz amostra×STR: `M[cbind(sample, STR)] <- 1`. NÃO usar
> `M[rows, cols]` (produto cartesiano).

## Como rodar (cluster, env `dbscan-r`)
```bash
qsub submit_burden.pbs        # Estratégia A (DBSCAN)
qsub submit_burden_gwas.pbs   # Estratégia B (GWAS-based)
```

## Saídas (`results/` e `results_gwas/`)
`burden_global.tsv`, `skat_per_gene.tsv`, `gene_burden.tsv`, `str_burden.tsv`,
`*_hits_uncorrected.tsv` / `*_hits_corrected.tsv` (BH q < 0.05).
