#!/usr/bin/env python3
"""
cross_DEGs_STRs.py
------------------
Cruza DEGs de multiplos datasets (GSE157103, GSE188847, GSE183533)
com o catalogo de STRs da coorte, incluindo metricas DBSCAN global e GWAS.

Saidas:
  1) all_STRs_in_DEGs.tsv       – todos os STRs anotados nos genes DEGs
  2) outlier_STRs_in_DEGs.tsv   – apenas STRs com outliers DBSCAN global
  3) rna_gene_strs.tsv          – pares gene x STR dos DEGs (analogo GWAS
                                  suggestive_gene_strs.tsv), com acesso GSE
  4) rna_outlier_genes.tsv      – STRs outlier por gene/DEG (analogo GWAS
                                  covid_suggestive_genes_with_outlier_STRs.tsv),
                                  com acesso GSE de origem
  5) rna_outlier_genes_by_study.tsv – contagem de STRs outlier por (gene, GSE)
   6) rna_summary_by_study.tsv        – resumo pos-estudo por (estudo, gene):
                                        STRs identificadas e se ha sobreposicao
                                        do maior alelo entre grupos
                                        (sim/nao/sem_dados)

Uso:
  python3 cross_DEGs_STRs.py \
      --deg-dir <diretorio_com_subpastas_GSE> \
      --str-catalog <caminho/STRs_analysis_dataset.tsv> \
      --gwas-outliers <caminho/suggestive_strs_outliers.tsv> \
      --out-dir <diretorio_de_saida>
"""
import argparse
import csv
import glob
import os
import sys


GENE_COL_CANDIDATES = ['Gene', 'gene', 'gene_name', 'gene_symbol']
SIGNIFICANT_COL_CANDIDATES = ['Significant', 'significance']
FDR_COL_CANDIDATES = ['FDR', 'adj.P.Val', 'P.Value']


def detect_gene_col(header):
    low = [h.strip().lower() for h in header]
    for cand in GENE_COL_CANDIDATES:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def detect_significant_col(header):
    low = [h.strip().lower() for h in header]
    for cand in SIGNIFICANT_COL_CANDIDATES:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def detect_fdr_col(header):
    low = [h.strip().lower() for h in header]
    for cand in FDR_COL_CANDIDATES:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def detect_logfc_col(header):
    low = [h.strip().lower() for h in header]
    for cand in ['logFC', 'logfc', 'log2fc']:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def detect_direction_col(header):
    low = [h.strip().lower() for h in header]
    for cand in ['Direction', 'direction', 'significance']:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def load_deg_file(path):
    genes = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        header = reader.fieldnames
        if not header:
            return genes

        gene_col = detect_gene_col(header)
        if gene_col is None:
            sys.stderr.write(f"  AVISO: coluna gene nao encontrada em {path}, pulando\n")
            return genes

        fdr_col = detect_fdr_col(header)
        logfc_col = detect_logfc_col(header)
        dir_col = detect_direction_col(header)
        sig_col = detect_significant_col(header)

        for row in reader:
            gene = row.get(gene_col, '').strip()
            if not gene:
                continue

            significant = True
            if sig_col:
                val = row.get(sig_col, '').strip().lower()
                significant = val in ('yes', 'true', '1', 'up', 'down',
                                      'up_covid', 'down_covid',
                                      'up_hfd0', 'down_hfd0',
                                      'up_icuvent', 'down_icuvent')

            if not significant:
                continue

            genes[gene] = {
                'logFC': row.get(logfc_col, '') if logfc_col else '',
                'FDR': row.get(fdr_col, '') if fdr_col else '',
                'Direction': row.get(dir_col, '') if dir_col else '',
            }
    return genes


def load_gwas_outliers(path):
    """Carrega suggestive_strs_outliers.tsv e retorna dict por STRs_ID."""
    gwas = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            sid = row.get('STRs_ID', '').strip()
            if not sid:
                continue
            gwas[sid] = {
                'n_outliers_dbscan_gwas': row.get('n_outliers', ''),
                'outlier_samples_dbscan_gwas': row.get('outlier_samples', ''),
                'outlier_residuals_dbscan_gwas': row.get('outlier_residuals', ''),
                'n_clusters_dbscan_gwas': row.get('n_clusters', ''),
                'noise_ratio_dbscan_gwas': row.get('noise_ratio', ''),
            }
    return gwas


def main():
    ap = argparse.ArgumentParser(
        description='Cruza DEGs de multiplos GSE com STRs da coorte')
    ap.add_argument('--deg-dir', required=True,
                    help='Diretorio raiz com subpastas GSE')
    ap.add_argument('--str-catalog', required=True,
                    help='Caminho para STRs_analysis_dataset.tsv')
    ap.add_argument('--gwas-outliers', required=True,
                    help='Caminho para suggestive_strs_outliers.tsv (DBSCAN GWAS)')
    ap.add_argument('--out-dir', default='.',
                    help='Diretorio de saida')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    gse_dirs = sorted([d for d in glob.glob(os.path.join(args.deg_dir, 'GSE*'))
                        if os.path.isdir(d)])
    if not gse_dirs:
        sys.stderr.write(f"ERRO: nenhuma pasta GSE encontrada em {args.deg_dir}\n")
        sys.exit(1)
    sys.stderr.write(f"Pastas GSE encontradas: {len(gse_dirs)}\n")

    sys.stderr.write(f"Carregando catalogo de STRs: {args.str_catalog}\n")
    str_catalog = []
    with open(args.str_catalog) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            str_catalog.append(row)
    sys.stderr.write(f"  {len(str_catalog)} STRs carregados\n")

    sys.stderr.write(f"Carregando DBSCAN GWAS: {args.gwas_outliers}\n")
    gwas_outliers = load_gwas_outliers(args.gwas_outliers)
    sys.stderr.write(f"  {len(gwas_outliers)} STRs com DBSCAN GWAS\n")

    all_matches = []
    outlier_matches = []

    for gse_dir in gse_dirs:
        gse_name = os.path.basename(gse_dir)
        deg_files = sorted(glob.glob(os.path.join(gse_dir, '*.tsv')) +
                           glob.glob(os.path.join(gse_dir, '*.csv')))

        for deg_file in deg_files:
            fname = os.path.basename(deg_file)
            dataset_label = f"{gse_name}/{fname}"
            sys.stderr.write(f"\n=== {dataset_label} ===\n")

            degs = load_deg_file(deg_file)
            sys.stderr.write(f"  {len(degs)} genes DEGs (Significant=Yes)\n")
            if not degs:
                continue

            n_all = 0
            n_outlier = 0
            for str_row in str_catalog:
                gene_name = str_row.get('gene_name', '').strip()
                if gene_name not in degs:
                    continue

                sid = str_row.get('STRs_ID', '').strip()
                gwas = gwas_outliers.get(sid, {})

                match = {
                    'dataset': dataset_label,
                    'gene_name': gene_name,
                    'STRs_ID': sid,
                    'chrom': str_row.get('chrom', ''),
                    'start': str_row.get('start', ''),
                    'end': str_row.get('end', ''),
                    'repeat_unit': str_row.get('repeat_unit', ''),
                    'allele1_est': str_row.get('allele1_est', ''),
                    'allele2_est': str_row.get('allele2_est', ''),
                    'depth': str_row.get('depth', ''),
                    'region': str_row.get('region', ''),
                    'group': str_row.get('group', ''),
                    'logFC': degs[gene_name]['logFC'],
                    'FDR': degs[gene_name]['FDR'],
                    'Direction': degs[gene_name]['Direction'],
                    'n_outliers_dbscan_global': str_row.get('n_outliers_dbscan_global', ''),
                    'outlier_samples_dbscan_global': str_row.get('outlier_samples_dbscan_global', ''),
                    'outlier_residuals_dbscan_global': str_row.get('outlier_residuals_dbscan_global', ''),
                    'n_clusters_dbscan_global': str_row.get('n_clusters_dbscan_global', ''),
                    'noise_ratio_dbscan_global': str_row.get('noise_ratio_dbscan_global', ''),
                    'n_outliers_dbscan_gwas': gwas.get('n_outliers_dbscan_gwas', ''),
                    'outlier_samples_dbscan_gwas': gwas.get('outlier_samples_dbscan_gwas', ''),
                    'outlier_residuals_dbscan_gwas': gwas.get('outlier_residuals_dbscan_gwas', ''),
                    'n_clusters_dbscan_gwas': gwas.get('n_clusters_dbscan_gwas', ''),
                    'noise_ratio_dbscan_gwas': gwas.get('noise_ratio_dbscan_gwas', ''),
                }
                all_matches.append(match)
                n_all += 1

                try:
                    n_clusters = int(match['n_clusters_dbscan_global'])
                    noise_ratio = float(match['noise_ratio_dbscan_global'])
                    n_outliers_val = int(match['n_outliers_dbscan_global'])
                except (ValueError, KeyError):
                    continue
                if n_clusters > 0 and noise_ratio <= 0.10 and n_outliers_val >= 1:
                    outlier_matches.append(match)
                    n_outlier += 1

            sys.stderr.write(f"  STRs encontrados: {n_all} total, "
                             f"{n_outlier} com outliers DBSCAN global\n")

    fields = ['dataset', 'gene_name', 'STRs_ID', 'chrom', 'start', 'end',
              'repeat_unit', 'allele1_est', 'allele2_est', 'depth', 'region', 'group',
              'logFC', 'FDR', 'Direction',
              'n_outliers_dbscan_global', 'outlier_samples_dbscan_global',
              'outlier_residuals_dbscan_global', 'n_clusters_dbscan_global',
              'noise_ratio_dbscan_global',
              'n_outliers_dbscan_gwas', 'outlier_samples_dbscan_gwas',
              'outlier_residuals_dbscan_gwas', 'n_clusters_dbscan_gwas',
              'noise_ratio_dbscan_gwas']

    out_all = os.path.join(args.out_dir, 'all_STRs_in_DEGs.tsv')
    with open(out_all, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(all_matches)
    sys.stderr.write(f"\nEscrito: {out_all} ({len(all_matches)} linhas)\n")

    out_outlier = os.path.join(args.out_dir, 'outlier_STRs_in_DEGs.tsv')
    with open(out_outlier, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(outlier_matches)
    sys.stderr.write(f"Escrito: {out_outlier} ({len(outlier_matches)} linhas)\n")

    # ----------------------------------------------------------------------
    # Saidas adicionais (análogas ao pipeline GWAS 6.4.1.2), com sufixo RNA
    # e indicando o estudo de origem (acesso GSE) de cada registro.
    # ----------------------------------------------------------------------

    def gse_name(dataset_label):
        # dataset = "GSE157103/<arquivo>..."; retorna o acesso GSE
        return dataset_label.split('/')[0]

    # 1) rna_gene_strs.tsv — background gene x STR dos DEGs (analogo de
    #    suggestive_gene_strs.tsv), 1 linha por (gene, strs_id), com datasets.
    pairs = {}   # (gene, STRs_ID) -> dict agregado
    for m in all_matches:
        key = (m['gene_name'], m['STRs_ID'])
        if key not in pairs:
            pairs[key] = {
                'gene': m['gene_name'],
                'strs_id': m['STRs_ID'],
                'chrom': m['chrom'],
                'start': m['start'],
                'end': m['end'],
                'repeat_unit': m['repeat_unit'],
                'region': m['region'],
                'datasets': [],
            }
        gse = gse_name(m['dataset'])
        if gse not in pairs[key]['datasets']:
            pairs[key]['datasets'].append(gse)

    bg_fields = ['gene', 'strs_id', 'chrom', 'start', 'end',
                 'repeat_unit', 'region', 'datasets']
    out_bg = os.path.join(args.out_dir, 'rna_gene_strs.tsv')
    with open(out_bg, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=bg_fields, delimiter='\t')
        w.writeheader()
        for key in sorted(pairs.keys()):
            row = dict(pairs[key])
            row['datasets'] = ';'.join(sorted(row['datasets']))
            w.writerow(row)
    sys.stderr.write(f"Escrito: {out_bg} ({len(pairs)} pares gene x STR)\n")

    # 2) rna_outlier_genes.tsv — análogo por gene de
    #    covid_suggestive_genes_with_outlier_STRs.tsv. 1 linha por
    #    (gene, strs_id, dataset) para DEGs com outlier DBSCAN global,
    #    incluindo a origem (acesso GSE).
    og_fields = ['gene', 'dataset', 'gse', 'strs_id', 'chrom', 'start', 'end',
                 'repeat_unit', 'region',
                 'logFC', 'FDR', 'Direction',
                 'n_outliers_dbscan_global', 'outlier_samples_dbscan_global',
                 'outlier_residuals_dbscan_global', 'n_clusters_dbscan_global',
                 'noise_ratio_dbscan_global',
                 'n_outliers_dbscan_gwas', 'outlier_samples_dbscan_gwas',
                 'outlier_residuals_dbscan_gwas', 'n_clusters_dbscan_gwas',
                 'noise_ratio_dbscan_gwas']
    out_og = os.path.join(args.out_dir, 'rna_outlier_genes.tsv')
    og_rows = []
    for m in outlier_matches:
        row = {
            'gene': m['gene_name'],
            'dataset': m['dataset'],
            'gse': gse_name(m['dataset']),
            'strs_id': m['STRs_ID'],
            'chrom': m['chrom'], 'start': m['start'], 'end': m['end'],
            'repeat_unit': m['repeat_unit'], 'region': m['region'],
            'logFC': m['logFC'], 'FDR': m['FDR'], 'Direction': m['Direction'],
            'n_outliers_dbscan_global': m['n_outliers_dbscan_global'],
            'outlier_samples_dbscan_global': m['outlier_samples_dbscan_global'],
            'outlier_residuals_dbscan_global': m['outlier_residuals_dbscan_global'],
            'n_clusters_dbscan_global': m['n_clusters_dbscan_global'],
            'noise_ratio_dbscan_global': m['noise_ratio_dbscan_global'],
            'n_outliers_dbscan_gwas': m['n_outliers_dbscan_gwas'],
            'outlier_samples_dbscan_gwas': m['outlier_samples_dbscan_gwas'],
            'outlier_residuals_dbscan_gwas': m['outlier_residuals_dbscan_gwas'],
            'n_clusters_dbscan_gwas': m['n_clusters_dbscan_gwas'],
            'noise_ratio_dbscan_gwas': m['noise_ratio_dbscan_gwas'],
        }
        og_rows.append(row)
    with open(out_og, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=og_fields, delimiter='\t')
        w.writeheader()
        w.writerows(og_rows)
    sys.stderr.write(f"Escrito: {out_og} ({len(og_rows)} linhas gene x STR x GSE)\n")

    # 3) rna_outlier_genes_per_datasets.tsv — resumo por (gene, GSE): quantos
    #    STRs outlier em cada estudo. Util para inspeção rapida.
    per_gse = {}
    for row in og_rows:
        key = (row['gene'], row['gse'])
        per_gse[key] = per_gse.get(key, 0) + 1
    out_pg = os.path.join(args.out_dir, 'rna_outlier_genes_by_study.tsv')
    with open(out_pg, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gene', 'gse', 'n_outlier_str_loci'])
        for (gene, gse), n in sorted(per_gse.items()):
            w.writerow([gene, gse, n])
    sys.stderr.write(f"Escrito: {out_pg} ({len(per_gse)} pares gene x GSE)\n")

    # 6) rna_summary_by_study.tsv — resumo pos-estudo por (estudo, gene):
    #    quantas STRs foram identificadas e se ha sobreposicao ou NAO do
    #    tamanho do alelo maior (maior valor entre allele1_est/allele2_est;
    #    na pratica allele2_est) entre grupos (coluna group do catalogo).
    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def largest_allele(m):
        a1 = num(m.get('allele1_est'))
        a2 = num(m.get('allele2_est'))
        if a1 is None and a2 is None:
            return None
        if a1 is None:
            return a2
        if a2 is None:
            return a1
        return max(a1, a2)

    # Valores do maior alelo por grupo, por STR locus (linhas = amostra x STR).
    locus_group_vals = {}
    for m in all_matches:
        sid = m['STRs_ID']
        g = (m.get('group') or '').strip()
        big = largest_allele(m)
        if not g or big is None:
            continue
        locus_group_vals.setdefault(sid, {}).setdefault(g, []).append(big)

    def locus_overlap(sid):
        groups = locus_group_vals.get(sid, {})
        labels = [g for g, vals in groups.items() if vals]
        if len(labels) < 2:
            return None
        ranges = {g: (min(groups[g]), max(groups[g])) for g in labels}
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a = ranges[labels[i]]
                b = ranges[labels[j]]
                if not (a[0] <= b[1] and b[0] <= a[1]):
                    return False
        return True

    # STR loci por gene (catálogo é da coorte, independe do estudo).
    gene_str_loci = {}
    for m in all_matches:
        gene_str_loci.setdefault(m['gene_name'], set()).add(m['STRs_ID'])

    def gene_overlap(gene):
        flags = [locus_overlap(sid) for sid in gene_str_loci.get(gene, ())]
        flags = [f for f in flags if f is not None]
        if not flags:
            return 'sem_dados'
        return 'nao' if any(f is False for f in flags) else 'sim'

    def gse_gene_str_loci(gse, rows):
        out = {}
        for m in rows:
            if gse_name(m['dataset']) != gse:
                continue
            out.setdefault(m['gene_name'], set()).add(m['STRs_ID'])
        return out

    all_gses = sorted({gse_name(m['dataset']) for m in all_matches}) if all_matches else []
    outlier_gene_pairs = {}
    for m in outlier_matches:
        outlier_gene_pairs.setdefault(gse_name(m['dataset']), set()) \
            .add((m['gene_name'], m['STRs_ID']))

    summary_by_study = []
    for gse in all_gses:
        gene_loci = gse_gene_str_loci(gse, all_matches)
        for gene in sorted(gene_loci):
            loci = gene_loci[gene]
            outlier_loci = {sid for (gn, sid) in outlier_gene_pairs.get(gse, ())
                            if gn == gene}
            summary_by_study.append({
                'gse': gse,
                'gene': gene,
                'n_strs_identified': len(loci),
                'n_strs_identified_outliers': len(outlier_loci),
                'overlap_maior_alealo_grupos': gene_overlap(gene),
            })

    out_sum = os.path.join(args.out_dir, 'rna_summary_by_study.tsv')
    with open(out_sum, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gse', 'gene', 'n_strs_identified',
                                          'n_strs_identified_outliers',
                                          'overlap_maior_alealo_grupos'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_by_study)
    sys.stderr.write(f"Escrito: {out_sum} ({len(summary_by_study)} linhas gene x estudo)\n")

    sys.stderr.write("\n=== Resumo por estudo (GSE) ===\n")
    for gse in all_gses:
        rows = [r for r in summary_by_study if r['gse'] == gse]
        n_strs = sum(r['n_strs_identified'] for r in rows)
        n_out = sum(r['n_strs_identified_outliers'] for r in rows)
        n_sem = sum(1 for r in rows if r['overlap_maior_alealo_grupos'] == 'nao')
        sys.stderr.write(
            f"  {gse}: {n_strs} STRs identificadas "
            f"({n_out} com outliers DBSCAN global), {len(rows)} genes; "
            f"{n_sem} gene(s) SEM sobreposicao do alelo maior\n")
    n_sem_nao = sum(1 for r in summary_by_study
                    if r['overlap_maior_alealo_grupos'] == 'nao')
    sys.stderr.write(
        f"  TOTAL: {len(summary_by_study)} pares gene x estudo; "
        f"{n_sem_nao} par(es) SEM sobreposicao do alelo maior\n")

    sys.stderr.write("\nConcluido.\n")


if __name__ == '__main__':
    main()
