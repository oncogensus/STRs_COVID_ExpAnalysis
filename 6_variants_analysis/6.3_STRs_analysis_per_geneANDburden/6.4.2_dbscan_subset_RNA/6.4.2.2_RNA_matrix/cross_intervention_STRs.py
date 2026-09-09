#!/usr/bin/env python3
"""
cross_intervention_STRs.py
--------------------------
Cruza DEGs de cada intervenção (arquivo TSV por comparação) com o catálogo
de STRs da coorte, incluindo métricas DBSCAN global.

Cada subpasta GSE contém um ou mais TSVs de DEGs. O nome da intervenção
é extraído do nome do arquivo (prefixo DEG(s)_ e sufixos removidos).

Saídas:
  1) intervention_strs.tsv        – todos os STRs anotados nos genes DEGs
  2) intervention_outliers.tsv    – apenas STRs com outliers DBSCAN global
  3) intervention_summary.tsv     – resumo por (intervenção, gene)

Uso:
  python3 cross_intervention_STRs.py \
      --deg-dir <diretorio_com_subpastas_GSE> \
      --str-catalog <caminho/STRs_analysis_dataset.tsv> \
      --out-dir <diretorio_de_saida>
"""
import argparse
import csv
import glob
import os
import re
import sys


GENE_COL_CANDIDATES = ['gene_symbol', 'Gene', 'gene', 'gene_name']
FDR_COL_CANDIDATES = ['FDR', 'adj.P.Val', 'P.Value']
LOGFC_COL_CANDIDATES = ['logFC', 'logfc', 'log2fc']
SIG_COL_CANDIDATES = ['Significant', 'significance']
DIR_COL_CANDIDATES = ['Direction', 'direction']


def detect_col(header, candidates):
    low = [h.strip().lower() for h in header]
    for cand in candidates:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def extract_intervention(fname):
    """Extrai nome da intervenção do nome do arquivo TSV."""
    base = os.path.splitext(fname)[0]
    # Remover prefixo GSEXXXXXX_ se existir
    base = re.sub(r'^GSE\d+_', '', base)
    # Remover prefixo DEG(s)_ 
    base = re.sub(r'^DEGs?_', '', base)
    # Remover sufixos conhecidos
    base = re.sub(r'_FDR0\.05_log2FC1$', '', base)
    base = re.sub(r'_ajustado$', '', base)
    base = re.sub(r'_FINAL$', '', base)
    base = re.sub(r'_FDR$', '', base)
    base = re.sub(r'_significativos$', '', base)
    base = re.sub(r'_significativo$', '', base)
    return base


def load_deg_file(path):
    genes = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        header = reader.fieldnames
        if not header:
            return genes

        gene_col = detect_col(header, GENE_COL_CANDIDATES)
        fdr_col = detect_col(header, FDR_COL_CANDIDATES)
        logfc_col = detect_col(header, LOGFC_COL_CANDIDATES)
        dir_col = detect_col(header, DIR_COL_CANDIDATES)
        sig_col = detect_col(header, SIG_COL_CANDIDATES)

        if gene_col is None:
            sys.stderr.write(f"  AVISO: coluna gene nao encontrada em {path}\n")
            return genes

        for row in reader:
            gene = row.get(gene_col, '').strip()
            if not gene:
                continue

            # Filtrar por significância
            if sig_col:
                val = row.get(sig_col, '').strip().lower()
                if val not in ('yes', 'true', '1', 'up', 'down',
                               'up_covid', 'down_covid',
                               'up_hfd0', 'down_hfd0',
                               'up_icuvent', 'down_icuvent'):
                    continue

            # Filtrar por FDR < 0.05 se disponível
            if fdr_col:
                try:
                    fdr_val = float(row.get(fdr_col, '1'))
                    if fdr_val >= 0.05:
                        continue
                except (ValueError, TypeError):
                    pass

            # Filtrar por |logFC| > 1 se disponível
            if logfc_col:
                try:
                    lfc = float(row.get(logfc_col, '0'))
                    if abs(lfc) < 1:
                        continue
                except (ValueError, TypeError):
                    pass

            genes[gene] = {
                'logFC': row.get(logfc_col, '') if logfc_col else '',
                'FDR': row.get(fdr_col, '') if fdr_col else '',
                'Direction': row.get(dir_col, '') if dir_col else '',
            }
    return genes


def main():
    ap = argparse.ArgumentParser(
        description='Cruza DEGs por intervenção com STRs da coorte')
    ap.add_argument('--deg-dir', required=True,
                    help='Diretório raiz com subpastas GSE')
    ap.add_argument('--str-catalog', required=True,
                    help='Caminho para STRs_analysis_dataset.tsv')
    ap.add_argument('--out-dir', default='.',
                    help='Diretório de saída')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Encontrar subpastas GSE
    gse_dirs = sorted([d for d in glob.glob(os.path.join(args.deg_dir, 'GSE*'))
                        if os.path.isdir(d)])
    if not gse_dirs:
        sys.stderr.write(f"ERRO: nenhuma pasta GSE encontrada em {args.deg_dir}\n")
        sys.exit(1)
    sys.stderr.write(f"Pastas GSE encontradas: {len(gse_dirs)}\n")

    # Carregar catálogo de STRs
    sys.stderr.write(f"Carregando catálogo de STRs: {args.str_catalog}\n")
    str_catalog = []
    with open(args.str_catalog) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            str_catalog.append(row)
    sys.stderr.write(f"  {len(str_catalog)} STRs carregados\n")

    all_matches = []
    outlier_matches = []
    intervention_info = {}  # intervention_name -> {gse, n_degs, n_files}

    for gse_dir in gse_dirs:
        gse_name = os.path.basename(gse_dir)
        deg_files = sorted(glob.glob(os.path.join(gse_dir, '*.tsv')))

        for deg_file in deg_files:
            fname = os.path.basename(deg_file)
            intervention = extract_intervention(fname)
            dataset_label = f"{gse_name}/{fname}"
            sys.stderr.write(f"\n=== {dataset_label} (intervenção: {intervention}) ===\n")

            degs = load_deg_file(deg_file)
            sys.stderr.write(f"  {len(degs)} genes DEGs\n")
            if not degs:
                continue

            intervention_info[intervention] = {
                'gse': gse_name,
                'file': fname,
                'n_degs': len(degs),
            }

            n_all = 0
            n_outlier = 0
            for str_row in str_catalog:
                gene_name = str_row.get('gene_name', '').strip()
                if gene_name not in degs:
                    continue

                sid = str_row.get('STRs_ID', '').strip()

                match = {
                    'intervention': intervention,
                    'gse': gse_name,
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
                    'n_clusters_dbscan_global': str_row.get('n_clusters_dbscan_global', ''),
                    'noise_ratio_dbscan_global': str_row.get('noise_ratio_dbscan_global', ''),
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

    # ----------------------------------------------------------------------
    # Saída 1: intervention_strs.tsv
    # ----------------------------------------------------------------------
    fields = ['intervention', 'gse', 'dataset', 'gene_name', 'STRs_ID',
              'chrom', 'start', 'end', 'repeat_unit',
              'allele1_est', 'allele2_est', 'depth', 'region', 'group',
              'logFC', 'FDR', 'Direction',
              'n_outliers_dbscan_global', 'outlier_samples_dbscan_global',
              'n_clusters_dbscan_global', 'noise_ratio_dbscan_global']

    out_strs = os.path.join(args.out_dir, 'intervention_strs.tsv')
    with open(out_strs, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(all_matches)
    sys.stderr.write(f"\nEscrito: {out_strs} ({len(all_matches)} linhas)\n")

    # ----------------------------------------------------------------------
    # Saída 2: intervention_outliers.tsv
    # ----------------------------------------------------------------------
    out_out = os.path.join(args.out_dir, 'intervention_outliers.tsv')
    with open(out_out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(outlier_matches)
    sys.stderr.write(f"Escrito: {out_out} ({len(outlier_matches)} linhas)\n")

    # ----------------------------------------------------------------------
    # Saída 3: intervention_summary.tsv (por intervenção x gene)
    # ----------------------------------------------------------------------
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

    # Valores do maior alelo por grupo, por STR locus
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

    # STR loci por gene (global, independe de intervenção)
    gene_str_loci = {}
    for m in all_matches:
        gene_str_loci.setdefault(m['gene_name'], set()).add(m['STRs_ID'])

    def gene_overlap(gene):
        flags = [locus_overlap(sid) for sid in gene_str_loci.get(gene, ())]
        flags = [f for f in flags if f is not None]
        if not flags:
            return 'sem_dados'
        return 'nao' if any(f is False for f in flags) else 'sim'

    # Agrupar por intervenção
    interventions = sorted(set(m['intervention'] for m in all_matches))

    # Pares gene x STR com outlier, por intervenção
    outlier_gene_pairs = {}
    for m in outlier_matches:
        outlier_gene_pairs.setdefault(m['intervention'], set()) \
            .add((m['gene_name'], m['STRs_ID']))

    summary_rows = []
    for interv in interventions:
        interv_matches = [m for m in all_matches if m['intervention'] == interv]
        gene_loci = {}
        for m in interv_matches:
            gene_loci.setdefault(m['gene_name'], set()).add(m['STRs_ID'])

        for gene in sorted(gene_loci):
            loci = gene_loci[gene]
            outlier_loci = {sid for (gn, sid) in outlier_gene_pairs.get(interv, ())
                            if gn == gene}
            summary_rows.append({
                'intervention': interv,
                'gse': intervention_info[interv]['gse'],
                'gene': gene,
                'n_strs_identified': len(loci),
                'n_strs_identified_outliers': len(outlier_loci),
                'overlap_maior_alealo_grupos': gene_overlap(gene),
            })

    out_sum = os.path.join(args.out_dir, 'intervention_summary.tsv')
    with open(out_sum, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['intervention', 'gse', 'gene',
                                          'n_strs_identified',
                                          'n_strs_identified_outliers',
                                          'overlap_maior_alealo_grupos'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)
    sys.stderr.write(f"Escrito: {out_sum} ({len(summary_rows)} linhas)\n")

    # Resumo
    sys.stderr.write("\n=== Resumo por intervenção ===\n")
    for interv in interventions:
        rows = [r for r in summary_rows if r['intervention'] == interv]
        n_strs = sum(r['n_strs_identified'] for r in rows)
        n_out = sum(r['n_strs_identified_outliers'] for r in rows)
        n_nao = sum(1 for r in rows if r['overlap_maior_alealo_grupos'] == 'nao')
        sys.stderr.write(
            f"  {interv}: {n_strs} STRs ({n_out} outliers), "
            f"{len(rows)} genes; {n_nao} SEM sobreposição\n")

    sys.stderr.write("\nConcluido.\n")


if __name__ == '__main__':
    main()
