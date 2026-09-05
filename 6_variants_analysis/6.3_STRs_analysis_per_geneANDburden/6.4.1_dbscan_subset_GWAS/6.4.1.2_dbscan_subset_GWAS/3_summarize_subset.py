#!/usr/bin/env python3
# summarize_subset.py
# ---------------------------------------------------------------------------
# PROPOSITO
#   Junta o dataset unificado (STRs_analysis_dataset_with_GWAS.tsv) com o
#   resultado do DBSCAN re-rodado no subset (suggestive_strs_outliers.tsv),
#   filtrando apenas os STRs com outlier no subset (n_outliers > 0).
#
#   Gera dois arquivos:
#     1) --out      Tabela final com anotacoes + outlier info
#     2) --summary  Resumo quantitativo: overlap, sinais DBSCAN, cobertura
#
# ENTRADAS
#   --unified    Dataset unificado (STRs_analysis_dataset_with_GWAS.tsv):
#                 todas as linhas do catalogo + gwas_hit, gwas_p, gwas_phenotypes,
#                 gwas_lead_snp, n_outliers_dbscan_global, etc.
#   --dbscan     Resultado do DBSCAN subset (suggestive_strs_outliers.tsv):
#                 STRs_ID, n_outliers, outlier_samples, outlier_residuals, etc.
#   --out        Tabela final anotada.
#   --summary    Resumo quantitativo (TSV).
# ---------------------------------------------------------------------------
import argparse, csv, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--unified', required=True)
    ap.add_argument('--dbscan', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--summary', required=True)
    args = ap.parse_args()

    # ----------------------------------------------------------------------
    # 1) Ler dataset unificado -> mapear info por STRs_ID
    #    Cada STRs_ID pode ter multiplas linhas (1 por sample).
    #    Interessa: gwas_hit, anotacoes, n_outliers_dbscan_global
    # ----------------------------------------------------------------------
    str_info = {}       # STRs_ID -> {gene, gwas_hit, gwas_p, ..., global_outliers}
    gwas_hit_genes = set()
    total_lines = 0
    total_strs = set()

    with open(args.unified) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            total_lines += 1
            sid = row.get('STRs_ID', '')
            if not sid:
                continue
            total_strs.add(sid)

            if sid not in str_info:
                # Pega anotacoes da primeira linha (gene_name, region, etc.)
                gwas_hit = row.get('gwas_hit', '0')
                try:
                    global_outliers = int(row.get('n_outliers_dbscan_global', '0') or '0')
                except ValueError:
                    global_outliers = 0

                str_info[sid] = {
                    'gene_name': row.get('gene_name', ''),
                    'region': row.get('region', ''),
                    'chrom': row.get('chrom', ''),
                    'start': row.get('start', ''),
                    'end': row.get('end', ''),
                    'repeat_unit': row.get('repeat_unit', ''),
                    'gwas_hit': gwas_hit,
                    'gwas_p': row.get('gwas_p', ''),
                    'gwas_phenotypes': row.get('gwas_phenotypes', ''),
                    'gwas_lead_snp': row.get('gwas_lead_snp', ''),
                    'global_outliers': global_outliers,
                }
                if gwas_hit == '1':
                    gwas_hit_genes.add(row.get('gene_name', ''))

    sys.stderr.write(f"Dataset unificado: {total_lines} linhas, "
                     f"{len(total_strs)} STRs_ID unicos\n")
    sys.stderr.write(f"Genes sugestivos com overlap: {len(gwas_hit_genes)}\n")

    # ----------------------------------------------------------------------
    # 2) Ler DBSCAN subset -> STRs com outlier no subset
    # ----------------------------------------------------------------------
    subset_outliers = {}  # STRs_ID -> (n_outliers, outlier_samples, ...)
    with open(args.dbscan) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            sid = row.get('STRs_ID', '')
            try:
                no = int(row.get('n_outliers', '0') or '0')
            except ValueError:
                no = 0
            if no > 0:
                subset_outliers[sid] = (no, row.get('outlier_samples', ''),
                                        row.get('outlier_residuals', ''),
                                        row.get('n_clusters', ''),
                                        row.get('noise_ratio', ''),
                                        row.get('eps', ''),
                                        row.get('minPts', ''),
                                        row.get('cutoff', ''),
                                        row.get('max_residual', ''),
                                        row.get('mean_residual', ''))

    sys.stderr.write(f"DBSCAN subset: {len(subset_outliers)} STRs com outlier (n_outliers >= 1)\n")

    # ----------------------------------------------------------------------
    # 3) Classificar STRs com GWAS hit por sinal DBSCAN
    # ----------------------------------------------------------------------
    # Contagens
    n_with_gwas_hit = 0
    n_with_global = 0
    n_with_subset = 0
    n_both = 0
    n_gwas_only = 0
    n_global_only = 0
    n_no_signal = 0

    # Para genes
    genes_with_global = set()
    genes_with_subset = set()
    genes_with_both = set()
    genes_with_gwas_only = set()
    genes_with_global_only = set()

    # Write output table
    out_rows = []
    for sid, info in str_info.items():
        if info['gwas_hit'] != '1':
            continue

        n_with_gwas_hit += 1
        gene = info['gene_name']
        has_global = info['global_outliers'] >= 1
        has_subset = sid in subset_outliers

        if has_global:
            n_with_global += 1
            genes_with_global.add(gene)
        if has_subset:
            n_with_subset += 1
            genes_with_subset.add(gene)

        if has_global and has_subset:
            n_both += 1
            genes_with_both.add(gene)
            signal = 'both'
        elif has_subset and not has_global:
            n_gwas_only += 1
            genes_with_gwas_only.add(gene)
            signal = 'gwas_only'
        elif has_global and not has_subset:
            n_global_only += 1
            genes_with_global_only.add(gene)
            signal = 'global_only'
        else:
            n_no_signal += 1
            signal = 'none'

        # Dados do subset DBSCAN (se existe)
        sub = subset_outliers.get(sid)
        sub_no = sub[0] if sub else 0
        sub_os = sub[1] if sub else ''
        sub_ore = sub[2] if sub else ''
        sub_nc = sub[3] if sub else ''
        sub_nr = sub[4] if sub else ''
        sub_eps = sub[5] if sub else ''
        sub_mpts = sub[6] if sub else ''
        sub_cut = sub[7] if sub else ''
        sub_mres = sub[8] if sub else ''
        sub_meres = sub[9] if sub else ''

        out_rows.append([
            gene, info['chrom'], info['start'], info['end'],
            info['gwas_p'], info['gwas_phenotypes'], info['gwas_lead_snp'],
            sid, info['region'], info['repeat_unit'],
            info['global_outliers'], signal,
            sub_no, sub_os, sub_ore, sub_nc, sub_nr,
            sub_eps, sub_mpts, sub_cut, sub_mres, sub_meres
        ])

    # Ordenar por gene, p-value
    out_rows.sort(key=lambda x: (x[0], float(x[4]) if x[4] else 1.0))

    # ----------------------------------------------------------------------
    # 4) Escrever tabela de saida
    # ----------------------------------------------------------------------
    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gene', 'chrom', 'gene_start', 'gene_end',
                     'gwas_p', 'gwas_phenotypes', 'gwas_lead_snp',
                     'strs_id', 'region', 'repeat_unit',
                     'global_outliers', 'signal_type',
                     'subset_n_outliers', 'subset_outlier_samples',
                     'subset_outlier_residuals', 'subset_n_clusters',
                     'subset_noise_ratio', 'subset_eps', 'subset_minPts',
                     'subset_cutoff', 'subset_max_residual', 'subset_mean_residual'])
        w.writerows(out_rows)

    sys.stderr.write(f"Escrevi {len(out_rows)} STRs com outlier no subset -> {args.out}\n")

    # ----------------------------------------------------------------------
    # 5) Gerar resumo quantitativo
    # ----------------------------------------------------------------------
    total_genes_suggestive = len(gwas_hit_genes)
    genes_with_any_signal = genes_with_both | genes_with_gwas_only | genes_with_global_only

    def pct(a, b):
        return f"{a/b*100:.1f}" if b > 0 else "0.0"

    summary_rows = [
        ['category', 'metric', 'value', 'percentage'],
        # --- Nivel STR ---
        ['STRs', 'total_strs_in_cohort', len(total_strs), '100.0'],
        ['STRs', 'strs_with_gwas_hit', n_with_gwas_hit,
         pct(n_with_gwas_hit, len(total_strs))],
        ['STRs', 'strs_with_global_dbscan_signal', n_with_global,
         pct(n_with_global, len(total_strs))],
        ['STRs', 'strs_with_subset_dbscan_signal', n_with_subset,
         pct(n_with_subset, len(total_strs))],
        ['STRs', 'strs_with_both_signals', n_both,
         pct(n_both, n_with_gwas_hit) if n_with_gwas_hit > 0 else '0.0'],
        ['STRs', 'strs_with_gwas_only_signal', n_gwas_only,
         pct(n_gwas_only, n_with_gwas_hit) if n_with_gwas_hit > 0 else '0.0'],
        ['STRs', 'strs_with_global_only_signal', n_global_only,
         pct(n_global_only, n_with_gwas_hit) if n_with_gwas_hit > 0 else '0.0'],
        ['STRs', 'strs_with_no_signal', n_no_signal,
         pct(n_no_signal, n_with_gwas_hit) if n_with_gwas_hit > 0 else '0.0'],
        # --- Nivel Gene ---
        ['Genes', 'total_suggestive_genes', total_genes_suggestive, '100.0'],
        ['Genes', 'genes_with_str_overlap', total_genes_suggestive,
         '100.0'],
        ['Genes', 'genes_with_any_dbscan_signal', len(genes_with_any_signal),
         pct(len(genes_with_any_signal), total_genes_suggestive)],
        ['Genes', 'genes_with_both_signals', len(genes_with_both),
         pct(len(genes_with_both), total_genes_suggestive)],
        ['Genes', 'genes_with_gwas_only', len(genes_with_gwas_only),
         pct(len(genes_with_gwas_only), total_genes_suggestive)],
        ['Genes', 'genes_with_global_only', len(genes_with_global_only),
         pct(len(genes_with_global_only), total_genes_suggestive)],
    ]

    with open(args.summary, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerows(summary_rows)

    sys.stderr.write(f"Resumo -> {args.summary}\n")

    # Imprime resumo no stderr
    sys.stderr.write("\n=== RESUMO ===\n")
    sys.stderr.write(f"STRs com overlap GWAS: {n_with_gwas_hit}\n")
    sys.stderr.write(f"  Sinal DBSCAN global:  {n_with_global} ({pct(n_with_global, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"  Sinal DBSCAN subset:  {n_with_subset} ({pct(n_with_subset, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"  Ambos:                {n_both} ({pct(n_both, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"  Apenas global:        {n_global_only} ({pct(n_global_only, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"  Apenas subset:        {n_gwas_only} ({pct(n_gwas_only, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"  Sem sinal:            {n_no_signal} ({pct(n_no_signal, n_with_gwas_hit)}%)\n")
    sys.stderr.write(f"\nGenes sugestivos: {total_genes_suggestive}\n")
    sys.stderr.write(f"  Com qualquer sinal DBSCAN: {len(genes_with_any_signal)} "
                     f"({pct(len(genes_with_any_signal), total_genes_suggestive)}%)\n")
    sys.stderr.write(f"  Com ambos:                {len(genes_with_both)} "
                     f"({pct(len(genes_with_both), total_genes_suggestive)}%)\n")
    sys.stderr.write(f"  Apenas global:            {len(genes_with_global_only)} "
                     f"({pct(len(genes_with_global_only), total_genes_suggestive)}%)\n")
    sys.stderr.write(f"  Apenas subset:            {len(genes_with_gwas_only)} "
                     f"({pct(len(genes_with_gwas_only), total_genes_suggestive)}%)\n")


if __name__ == '__main__':
    main()
