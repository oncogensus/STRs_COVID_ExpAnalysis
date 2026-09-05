#!/usr/bin/env python3
# summarize_subset.py
# ---------------------------------------------------------------------------
# PROPOSITO
#   Cruza o dataset unificado (STRs_analysis_dataset_with_GWAS.tsv) com o
#   resultado do DBSCAN subset (suggestive_strs_outliers.tsv) e gera um
#   resumo comparativo entre genes sugestivos (p<1e-5) e significativos (p<5e-8).
#
# SAIDAS
#   --out       Tabela final com anotacoes + sinal DBSCAN
#   --summary   Resumo quantitativo comparativo (sugestivo vs significativo)
#
# ENTRADAS
#   --unified           Dataset unificado (STRs_analysis_dataset_with_GWAS.tsv)
#   --dbscan            DBSCAN subset (suggestive_strs_outliers.tsv)
#   --covid-genes-sig   Genes significativos (covid_genes_significant.tsv, p<5e-8)
#   --out               Tabela final anotada
#   --summary           Resumo comparativo (TSV)
# ---------------------------------------------------------------------------
import argparse, csv, sys


def pct(a, b):
    return f"{a/b*100:.1f}" if b > 0 else "0.0"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--unified', required=True)
    ap.add_argument('--dbscan', required=True)
    ap.add_argument('--covid-genes-sig', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--summary', required=True)
    args = ap.parse_args()

    # ----------------------------------------------------------------------
    # 1) Ler genes significativos (p<5e-8) -> set de nomes
    # ----------------------------------------------------------------------
    sig_genes = set()
    with open(args.covid_genes_sig) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            gene = row.get('gene', '').strip()
            if gene:
                sig_genes.add(gene)
    sys.stderr.write(f"Genes significativos (p<5e-8): {len(sig_genes)}\n")

    # ----------------------------------------------------------------------
    # 2) Ler dataset unificado -> info por STRs_ID
    # ----------------------------------------------------------------------
    str_info = {}
    total_strs = set()
    gwas_hit_genes = set()
    n_total_rows = 0

    with open(args.unified) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            sid = row.get('STRs_ID', '')
            if not sid:
                continue
            n_total_rows += 1
            total_strs.add(sid)
            if sid not in str_info:
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

    # ----------------------------------------------------------------------
    # 3) Ler DBSCAN subset -> STRs com outlier
    # ----------------------------------------------------------------------
    subset_outliers = {}
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

    # ----------------------------------------------------------------------
    # 4) Contar outliers globais (todos os STRs, nao so com GWAS hit)
    # ----------------------------------------------------------------------
    total_strs_with_global_outliers = 0
    total_global_outlier_instances = 0
    strs_with_global_in_suggestive = 0
    global_outliers_in_suggestive = 0
    strs_with_global_in_significant = 0
    global_outliers_in_significant = 0

    for sid, info in str_info.items():
        n_out = info['global_outliers']
        if n_out >= 1:
            total_strs_with_global_outliers += 1
            total_global_outlier_instances += n_out
            if info['gwas_hit'] == '1':
                strs_with_global_in_suggestive += 1
                global_outliers_in_suggestive += n_out
                if info['gene_name'] in sig_genes:
                    strs_with_global_in_significant += 1
                    global_outliers_in_significant += n_out

    # ----------------------------------------------------------------------
    # 5) Classificar STRs e escrever tabela de saida
    # ----------------------------------------------------------------------
    # Contagens para sugestivo
    n_with_gwas_hit = 0
    n_with_global = 0
    n_with_subset = 0
    n_both = 0
    n_gwas_only = 0
    n_global_only = 0
    n_no_signal = 0

    # Contagens para significativo
    n_sig_with_hit = 0
    n_sig_with_global = 0
    n_sig_with_subset = 0
    n_sig_both = 0
    n_sig_gwas_only = 0
    n_sig_global_only = 0
    n_sig_no_signal = 0

    # Sets de genes
    genes_with_any_signal = set()
    genes_with_both = set()
    genes_with_gwas_only = set()
    genes_with_global_only = set()
    sig_genes_with_any_signal = set()
    sig_genes_with_both = set()
    sig_genes_with_gwas_only = set()
    sig_genes_with_global_only = set()

    out_rows = []
    for sid, info in str_info.items():
        if info['gwas_hit'] != '1':
            continue

        n_with_gwas_hit += 1
        gene = info['gene_name']
        is_sig = gene in sig_genes
        has_global = info['global_outliers'] >= 1
        has_subset = sid in subset_outliers

        if has_global:
            n_with_global += 1
            genes_with_global.add(gene) if 'genes_with_global' in dir() else None
        if has_subset:
            n_with_subset += 1

        if is_sig:
            n_sig_with_hit += 1
            if has_global:
                n_sig_with_global += 1
            if has_subset:
                n_sig_with_subset += 1

        if has_global and has_subset:
            n_both += 1
            genes_with_both.add(gene)
            signal = 'both'
            if is_sig:
                n_sig_both += 1
                sig_genes_with_both.add(gene)
        elif has_subset and not has_global:
            n_gwas_only += 1
            genes_with_gwas_only.add(gene)
            signal = 'subset_only'
            if is_sig:
                n_sig_gwas_only += 1
                sig_genes_with_gwas_only.add(gene)
        elif has_global and not has_subset:
            n_global_only += 1
            genes_with_global_only.add(gene)
            signal = 'global_only'
            if is_sig:
                n_sig_global_only += 1
                sig_genes_with_global_only.add(gene)
        else:
            n_no_signal += 1
            signal = 'none'

        genes_with_any_signal.add(gene)
        if is_sig:
            sig_genes_with_any_signal.add(gene)

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

        gwas_significance = 'significant' if is_sig else 'suggestive_only'

        out_rows.append([
            gene, info['chrom'], info['start'], info['end'],
            info['gwas_p'], gwas_significance, info['gwas_phenotypes'], info['gwas_lead_snp'],
            sid, info['region'], info['repeat_unit'],
            info['global_outliers'], signal,
            sub_no, sub_os, sub_ore, sub_nc, sub_nr,
            sub_eps, sub_mpts, sub_cut, sub_mres, sub_meres
        ])

    out_rows.sort(key=lambda x: (x[0], float(x[4]) if x[4] else 1.0))

    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gene', 'chrom', 'gene_start', 'gene_end',
                     'gwas_p', 'gwas_significance', 'gwas_phenotypes', 'gwas_lead_snp',
                     'strs_id', 'region', 'repeat_unit',
                     'global_outliers', 'signal_type',
                     'subset_n_outliers', 'subset_outlier_samples',
                     'subset_outlier_residuals', 'subset_n_clusters',
                     'subset_noise_ratio', 'subset_eps', 'subset_minPts',
                     'subset_cutoff', 'subset_max_residual', 'subset_mean_residual'])
        w.writerows(out_rows)

    sys.stderr.write(f"Escrevi {len(out_rows)} STRs com outlier no subset -> {args.out}\n")

    # ----------------------------------------------------------------------
    # 5) Quebra por fenotipo (A2, B2, C2 e combinacoes)
    # ----------------------------------------------------------------------
    from collections import defaultdict
    pheno_genes = defaultdict(set)       # pheno -> set(genes)
    pheno_strs = defaultdict(set)        # pheno -> set(strs_id)
    pheno_global = defaultdict(set)      # pheno -> set(strs_id com global outlier)
    pheno_subset = defaultdict(set)      # pheno -> set(strs_id com subset outlier)
    pheno_genes_with_signal = defaultdict(set)  # pheno -> set(genes com qualquer sinal DBSCAN)

    for sid, info in str_info.items():
        if info['gwas_hit'] != '1':
            continue
        phenos_raw = info.get('gwas_phenotypes', '')
        if not phenos_raw:
            continue
        # Fenotipos podem ser "A2" ou "A2,B2" ou "A2,B2,C2"
        phenos = [p.strip() for p in phenos_raw.split(',') if p.strip()]
        has_global = info['global_outliers'] >= 1
        has_subset = sid in subset_outliers
        has_signal = has_global or has_subset
        for ph in phenos:
            pheno_genes[ph].add(info['gene_name'])
            pheno_strs[ph].add(sid)
            if has_global:
                pheno_global[ph].add(sid)
            if has_subset:
                pheno_subset[ph].add(sid)
            if has_signal:
                pheno_genes_with_signal[ph].add(info['gene_name'])

    # Ordena fenotipos por nome
    sorted_phenos = sorted(pheno_genes.keys())

    # ----------------------------------------------------------------------
    # 6) Resumo comparativo
    # ----------------------------------------------------------------------
    total_genes_sugg = len(gwas_hit_genes)
    total_genes_sig = len(sig_genes)

    w = None
    with open(args.summary, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['section', 'metric', 'suggestive', 'suggestive_pct',
                     'significant', 'significant_pct', 'denominator'])

        def add(section, metric, sugg_val, sig_val, denom_sugg, denom_sig):
            w.writerow([section, metric, sugg_val,
                        pct(sugg_val, denom_sugg),
                        sig_val,
                        pct(sig_val, denom_sig),
                        f"sugg={denom_sugg};sig={denom_sig}"])

        add('genes', 'total', total_genes_sugg, total_genes_sig,
            total_genes_sugg, total_genes_sig)
        add('genes', 'with_str_in_cohort', total_genes_sugg, total_genes_sig,
            total_genes_sugg, total_genes_sig)
        add('overlap', 'strs_with_gwas_hit', n_with_gwas_hit, n_sig_with_hit,
            len(total_strs), len(total_strs))
        add('signal_strs', 'global_only', n_global_only, n_sig_global_only,
            n_with_gwas_hit, n_sig_with_hit)
        add('signal_strs', 'subset_only', n_gwas_only, n_sig_gwas_only,
            n_with_gwas_hit, n_sig_with_hit)
        add('signal_strs', 'both', n_both, n_sig_both,
            n_with_gwas_hit, n_sig_with_hit)
        add('signal_strs', 'none', n_no_signal, n_sig_no_signal,
            n_with_gwas_hit, n_sig_with_hit)
        add('signal_genes', 'any_signal', len(genes_with_any_signal),
            len(sig_genes_with_any_signal), total_genes_sugg, total_genes_sig)
        add('signal_genes', 'both', len(genes_with_both),
            len(sig_genes_with_both), total_genes_sugg, total_genes_sig)
        add('signal_genes', 'global_only', len(genes_with_global_only),
            len(sig_genes_with_global_only), total_genes_sugg, total_genes_sig)
        add('signal_genes', 'subset_only', len(genes_with_gwas_only),
            len(sig_genes_with_gwas_only), total_genes_sugg, total_genes_sig)
        # Outliers globais
        add('global_outliers', 'strs_with_outliers', total_strs_with_global_outliers,
            total_strs_with_global_outliers, len(total_strs), len(total_strs))
        add('global_outliers', 'total_instances', total_global_outlier_instances,
            total_global_outlier_instances, len(total_strs), len(total_strs))
        add('global_outliers', 'in_suggestive_genes', strs_with_global_in_suggestive,
            strs_with_global_in_significant, total_strs_with_global_outliers,
            total_strs_with_global_outliers)
        add('global_outliers', 'instances_in_suggestive', global_outliers_in_suggestive,
            global_outliers_in_significant, total_global_outlier_instances,
            total_global_outlier_instances)
        # Fenotipos
        for ph in sorted_phenos:
            w.writerow(['phenotype', ph,
                         len(pheno_genes[ph]), '',
                         '', '',
                         f"genes={len(pheno_genes[ph])};strs={len(pheno_strs[ph])};"
                         f"global={len(pheno_global[ph])};subset={len(pheno_subset[ph])}"])

    # Imprime resumo formatado no stderr
    sep = "=" * 55
    sys.stderr.write(f"\n{sep}\n")
    sys.stderr.write(" RESUMO: Overlap GWAS x DBSCAN\n")
    sys.stderr.write(f"{sep}\n\n")

    sys.stderr.write("[Dataset]\n")
    sys.stderr.write(f"  Total linhas (sample × STR):     {n_total_rows:>10,}\n")
    sys.stderr.write(f"  STRs_ID unicos:                  {len(total_strs):>10,}\n\n")

    sys.stderr.write("[Genes COVID-19 HG]\n")
    sys.stderr.write(f"  {'':30s} {'Sugestivo':>12s} {'Significativo':>14s}\n")
    sys.stderr.write(f"  {'Total genes':30s} {total_genes_sugg:>12d} {total_genes_sig:>14d}\n")
    sys.stderr.write(f"  {'Com STR na coorte':30s} {total_genes_sugg:>12d} {total_genes_sig:>14d}\n\n")

    sys.stderr.write("[Overlap STR x Gene]\n")
    sys.stderr.write(f"  STRs com GWAS hit:         {n_with_gwas_hit:>8d}  (signif: {n_sig_with_hit})\n\n")

    sys.stderr.write("[Sinal DBSCAN — entre STRs com GWAS hit]\n")
    sys.stderr.write(f"  {'':30s} {'Sugestivo':>12s} {'Significativo':>14s}\n")
    sys.stderr.write(f"  {'Com sinal global':30s} {n_with_global:>12d} {n_sig_with_global:>14d}\n")
    sys.stderr.write(f"  {'Com sinal subset':30s} {n_with_subset:>12d} {n_sig_with_subset:>14d}\n")
    sys.stderr.write(f"  {'Ambos':30s} {n_both:>12d} {n_sig_both:>14d}\n")
    sys.stderr.write(f"  {'Apenas global':30s} {n_global_only:>12d} {n_sig_global_only:>14d}\n")
    sys.stderr.write(f"  {'Apenas subset':30s} {n_gwas_only:>12d} {n_sig_gwas_only:>14d}\n")
    sys.stderr.write(f"  {'Sem sinal':30s} {n_no_signal:>12d} {n_sig_no_signal:>14d}\n\n")

    sys.stderr.write("[Genes com sinal DBSCAN]\n")
    sys.stderr.write(f"  {'':30s} {'Sugestivo':>12s} {'Significativo':>14s}\n")
    sys.stderr.write(f"  {'Com qualquer sinal':30s} {len(genes_with_any_signal):>12d} {len(sig_genes_with_any_signal):>14d}\n")
    sys.stderr.write(f"  {'  % do total':30s} {pct(len(genes_with_any_signal), total_genes_sugg):>11s}% {pct(len(sig_genes_with_any_signal), total_genes_sig):>13s}%\n")
    sys.stderr.write(f"  {'Ambos':30s} {len(genes_with_both):>12d} {len(sig_genes_with_both):>14d}\n")
    sys.stderr.write(f"  {'Apenas global':30s} {len(genes_with_global_only):>12d} {len(sig_genes_with_global_only):>14d}\n")
    sys.stderr.write(f"  {'Apenas subset':30s} {len(genes_with_gwas_only):>12d} {len(sig_genes_with_gwas_only):>14d}\n\n")

    sys.stderr.write("[Outliers Globais — DBSCAN global (todos os STRs)]\n")
    sys.stderr.write(f"  Total STRs no cohort:                         {len(total_strs):>10,}\n")
    sys.stderr.write(f"  STRs com outlier global:                         {total_strs_with_global_outliers:>6,}  ({pct(total_strs_with_global_outliers, len(total_strs))}% do total)\n")
    sys.stderr.write(f"  Total instâncias outlier (sample × STR):         {total_global_outlier_instances:>6,}\n\n")
    sys.stderr.write(f"  Distribuição por gene:\n")
    sys.stderr.write(f"    Em genes sugestivos (p<1e-5):     {strs_with_global_in_suggestive:>4d} STRs  ({global_outliers_in_suggestive:>4d} instâncias)  [{pct(strs_with_global_in_suggestive, total_strs_with_global_outliers)}% dos outliers]\n")
    sys.stderr.write(f"    Em genes significativos (p<5e-8):  {strs_with_global_in_significant:>4d} STRs  ({global_outliers_in_significant:>4d} instâncias)  [{pct(strs_with_global_in_significant, total_strs_with_global_outliers)}% dos outliers]\n")
    n_outside = total_strs_with_global_outliers - strs_with_global_in_suggestive
    n_outside_inst = total_global_outlier_instances - global_outliers_in_suggestive
    sys.stderr.write(f"    Fora de genes sugestivos:         {n_outside:>4d} STRs  ({n_outside_inst:>4d} instâncias)  [{pct(n_outside, total_strs_with_global_outliers)}% dos outliers]\n\n")

    sys.stderr.write("[Fenótipos — genes sugestivos associados a cada fenótipo]\n")
    sys.stderr.write(f"  Cada linha contém apenas genes cujo gwas_phenotypes inclui aquele fenótipo.\n\n")
    sys.stderr.write(f"  {'Fenótipo':8s} {'Genes':>6s} {'STRs':>6s} {'Global':>7s} {'Subset':>7s}  {'Genes com sinal DBSCAN'}\n")
    sys.stderr.write(f"  {'─'*8} {'─'*6} {'─'*6} {'─'*7} {'─'*7}  {'─'*30}\n")
    for ph in sorted_phenos:
        gene_list = sorted(pheno_genes_with_signal.get(ph, set()))
        if len(gene_list) > 5:
            gene_str = ", ".join(gene_list[:5]) + f", ... (+{len(gene_list)-5})"
        elif gene_list:
            gene_str = ", ".join(gene_list)
        else:
            gene_str = "nenhum"
        sys.stderr.write(f"  {ph:8s} {len(pheno_genes[ph]):>6d} {len(pheno_strs[ph]):>6d} "
                         f"{len(pheno_global[ph]):>7d} {len(pheno_subset[ph]):>7d}  {gene_str}\n")

    sys.stderr.write(f"\n{sep}\n")
    sys.stderr.write(f"  Saida: {args.out}\n")
    sys.stderr.write(f"  Resumo: {args.summary}\n")
    sys.stderr.write(f"{sep}\n")


if __name__ == '__main__':
    main()
