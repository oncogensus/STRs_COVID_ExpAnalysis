# 6.4.7 — SKAT-O por gene/via + enriquecimento de vias (KEGG/Reactome)

Testa se STRs associam com COVID-19 (caso x controle) usando **SKAT-O**
(combinação ótima de teste *burden* + SKAT, com `rho` estimado), tanto por
**gene** quanto agrupados por **via KEGG/Reactome**, e depois mede **enriquecimento
de vias** para os genes candidatos de cada método (SKAT-O, burden, outliers
GWAS/RNA) e para os genes de **expressão diferencial (RNA-seq)**.

## Pipeline

| Passo | Script | O que faz | Saídas |
|-------|--------|-----------|--------|
| 0 | `00_common.R` | Funções compartilhadas (inputs, matriz de dosagem, modelo nulo, SKAT.O, mapas de vias) | — |
| 1 | `1_skat_o_per_gene.R` | SKAT.O por gene para `gwas_burden` e `rna_burden` | `results_skat_o_genes{, _rna}/skat_o_per_gene*.tsv`, `burden_per_gene*.tsv` |
| 2 | `2_skat_o_per_pathway.R` | SKAT.O por via KEGG e Reactome (mesmas estratégias) | `results_skat_o_pathways/skat_o_pathways_{kegg,reactome,all}.tsv`, `pathway_gene_*.tsv` |
| 3 | `3_enrichment_clusterprofiler.R` (+ `.pbs`, `3_run_enrichment_for_sources.sh`) | Enriquecimento de vias nos genes de cada fonte (SKAT-O, burden, outliers) via **clusterProfiler** (`enrichKEGG`/`enrichPathway`), com fallback hipergeométrico (msigdbr); cruza com via KEGG hsa05171 | `results_enrichment/enrichment_<label>.tsv`, `dotplot_<label>.png`, `covid19_pathway_overlap_<label>.tsv` |
| 4 | `4_de_enrichment.R` (+ `.pbs`) | Fisher DE × STR por estudo RNA-seq + vias enriquecidas em genes DE | `results_de_enrichment/de_str_fisher_by_study.tsv`, `de_pathway_enrichment.tsv` |

Rodar tudo: `qsub submit_all.sh` (4 jobs PBS) ou `bash run_all.sh` (login node).

## Modelo estatístico

- `SKAT_Null_Model(group ~ age + sex + EV1 + EV2 + EV3, out_type="D")` — caso=1/controle=0,
  ajuste por idade, sexo e 3 PCs de ancestralidade (mesmo do `burden_gwas.R` da 6.4.3).

### Comando do teste e compatibilidade de API do pacote SKAT

- `un_teste` chama `run_skat_o()` (do `00_common.R`), que roda o teste combinado
  SKAT-O sobre a matriz de STRs (`kernel="linear.weighted"`) e relata: `p_value_skat_o`
  (combinado), `p_value_skat`, `p_value_burden` e `rho` (peso ótimo entre os dois).
- **Versões do SKAT** (o pacote mudou a API nas versões ≥ 2.2):
  - SKAT ≤ 2.1: chamada direta `SKAT.O(Z, obj, kernel)` → `$p.value`, `$p.value.SKAT`,
    `$p.value.burden`, `$rho`.
  - SKAT ≥ 2.2: a função `SKAT.O` foi removida; o código cai automaticamente para
    `SKAT(Z, obj, kernel, r.corr = seq(0, 1, by = 0.1))`, extraindo o p combinado de
    `$p.value`, os p por `rho` de `$param$p.val.each` (1ª = SKAT, rho=0; última = burden,
    rho=1) e o `rho` ótimo de `$param$rho_est`. A coluna `note` registra qual API foi usada.
- Correção múltipla **BH** dentro de cada teste (gene, via-KEGG, via-Reactome).

## Estratégias de background

| Estratégia | STRs testados |
|------------|---------------|
| `gwas_burden` | STRs de genes GWAS-sugestivos (`6.4.1.2.../results/suggestive_gene_strs.tsv`) |
| `rna_burden`  | STRs de genes DEGs (`6.4.2.2.../results/rna_gene_strs.tsv`) |
| `full`        | todos os STRs com dosagem no `norm_file` |

## Mapas de vias

- **Passo 2**: obtém anotações via **clusterProfiler** (`download_KEGG("hsa")`,
  `download_Reactome("hsa")`), convertendo os IDs Entrez→Ensembl com `org.Hs.eg.db`
  (default), ou via **msigdbr** (C2 CP:KEGG / CP:REACTOME) com `--gene-set-source msigdbr`.
  O mapa via×gene é **cacheado em TSV** (`--cache-dir`, default `results_skat_o_pathways/pathway_cache`),
  então vias não precisam ser re-baixadas a cada rodada.
- **Passo 3**: preferencialmente `clusterProfiler::enrichKEGG`/`ReactomePA::enrichPathway`
  (Entrez via `bitr`); se esses pacotes não estiverem disponíveis (ou sem internet no
  nodo), cai para **teste hipergeométrico (`phyper`)** usando os mapas msigdbr cacheados
  no mesmo espaço Ensembl, com BH.

## Fontes de genes para enriquecimento (passo 3)

Extraídas automaticamente pelo `3_run_enrichment_for_sources.sh`:
- `skato_gwas_burden`, `skato_rna_burden` — genes com SKAT-O **q<0.05** (passo 1)
- `burden_gwas_unc`, `burden_rna_unc` — genes com burden **p<0.05** (passo 1)
- `outlier_union` — genes com STR outlier GWAS/RNA do comparativo 6.4.4

## Dependências (env `r_enrich_env`)

`Rscript` com: `SKAT`, `data.table`, `dplyr`, `ggplot2`, `clusterProfiler`,
`ReactomePA`, `org.Hs.eg.db`, `AnnotationDbi`, `msigdbr`.

## Inputs (paths fixos, como na 6.4.3)

- `5_global_dbscan/norm_test/STRs_normalized_residuals.tsv`
- `samples/samples_infos.csv`, `samples/STRs_analysis_dataset.tsv`
- `4_ancestry/EthSEQ_Results_3D/Report.PCAcoord`
- backgrounds `suggestive_gene_strs.tsv` / `rna_gene_strs.tsv`
- DEGs por estudo: edite `DEG_DIR` em `4_de_enrichment.pbs`

## Próximos ajustes sugeridos

- Definir universo explícito (`--universe-file`) no passo 3 se quiser enriquecimento
  contra o conjunto exato de genes testados (default: genes de interesse).
- No passo 2, se `full` for muito pesado, filtrar por `--min-strs` (default 2 por via).