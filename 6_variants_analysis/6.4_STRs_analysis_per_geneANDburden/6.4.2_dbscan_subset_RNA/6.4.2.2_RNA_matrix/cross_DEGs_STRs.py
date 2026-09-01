#!/usr/bin/env python3
"""
cross_DEGs_STRs.py
------------------
Cruza DEGs de multiplos datasets (GSE157103, GSE188847, GSE183533)
com o catalogo de STRs da coorte.

Saidas:
  1) all_STRs_in_DEGs.tsv       – todos os STRs anotados nos genes DEGs
  2) outlier_STRs_in_DEGs.tsv   – apenas STRs com outliers DBSCAN

Uso:
  python3 cross_DEGs_STRs.py \
      --deg-dir <diretorio_com_subpastas_GSE> \
      --str-catalog <caminho/STRs_analysis_dataset.tsv> \
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


def main():
    ap = argparse.ArgumentParser(
        description='Cruza DEGs de multiplos GSE com STRs da coorte')
    ap.add_argument('--deg-dir', required=True,
                    help='Diretorio raiz com subpastas GSE')
    ap.add_argument('--str-catalog', required=True,
                    help='Caminho para STRs_analysis_dataset.tsv')
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

                match = {
                    'dataset': dataset_label,
                    'gene_name': gene_name,
                    'STRs_ID': str_row.get('STRs_ID', ''),
                    'chrom': str_row.get('chrom', ''),
                    'start': str_row.get('start', ''),
                    'end': str_row.get('end', ''),
                    'repeat_unit': str_row.get('repeat_unit', ''),
                    'allele1_est': str_row.get('allele1_est', ''),
                    'allele2_est': str_row.get('allele2_est', ''),
                    'depth': str_row.get('depth', ''),
                    'region': str_row.get('region', ''),
                    'n_outliers': str_row.get('n_outliers', '0'),
                    'n_clusters': str_row.get('n_clusters', '0'),
                    'noise_ratio': str_row.get('noise_ratio', ''),
                    'outlier_samples': str_row.get('outlier_samples', ''),
                    'outlier_residuals': str_row.get('outlier_residuals', ''),
                    'group': str_row.get('group', ''),
                    'logFC': degs[gene_name]['logFC'],
                    'FDR': degs[gene_name]['FDR'],
                    'Direction': degs[gene_name]['Direction'],
                }
                all_matches.append(match)
                n_all += 1

                try:
                    n_clusters = int(match['n_clusters'])
                    noise_ratio = float(match['noise_ratio'])
                    n_outliers_val = int(match['n_outliers'])
                except (ValueError, KeyError):
                    continue
                if n_clusters > 0 and noise_ratio <= 0.10 and n_outliers_val >= 1:
                    outlier_matches.append(match)
                    n_outlier += 1

            sys.stderr.write(f"  STRs encontrados: {n_all} total, "
                             f"{n_outlier} com outliers DBSCAN\n")

    fields = ['dataset', 'gene_name', 'STRs_ID', 'chrom', 'start', 'end',
              'repeat_unit', 'allele1_est', 'allele2_est', 'depth', 'region',
              'n_outliers', 'n_clusters', 'noise_ratio',
              'outlier_samples', 'outlier_residuals', 'group',
              'logFC', 'FDR', 'Direction']

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

    sys.stderr.write("\n=== Resumo por gene ===\n")
    genes_summary = {}
    for m in all_matches:
        g = m['gene_name']
        if g not in genes_summary:
            genes_summary[g] = {'total': 0, 'outliers': 0, 'datasets': set()}
        genes_summary[g]['total'] += 1
        genes_summary[g]['datasets'].add(m['dataset'])
    for m in outlier_matches:
        g = m['gene_name']
        if g in genes_summary:
            genes_summary[g]['outliers'] += 1

    for gene in sorted(genes_summary.keys()):
        info = genes_summary[gene]
        ds = ', '.join(sorted(info['datasets']))
        sys.stderr.write(f"  {gene}: {info['total']} STRs, "
                         f"{info['outliers']} outliers ({ds})\n")

    sys.stderr.write("\nConcluido.\n")


if __name__ == '__main__':
    main()
