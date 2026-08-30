#!/usr/bin/env python3
# get_outlier_genes.py
# Do catálogo completo de STRs da coorte (STRs_analysis_dataset.tsv), seleciona os
# loci (STRs_ID) com n_outliers > 0 (marcados como outliers pelo DBSCAN global em
# 5_dbscan) e extrai o conjunto de genes que os contêm (gene_name).
# Saída:
#   data/outlier_genes.txt  (gene_name, n_outlier_str_loci)  -- uma linha por gene
#
# Uso:
#   python3 get_outlier_genes.py --catalog $CATALOG --out data/outlier_genes.txt
import argparse, csv, sys


def detect_cols(header):
    low = [h.lower().lstrip('#') for h in header]
    idx = {}
    for i, h in enumerate(low):
        if 'strs_id' in h and 'strs_id' not in idx:
            idx['strs_id'] = i
        if h in ('gene_name', 'gene') and 'gene_name' not in idx:
            idx['gene_name'] = i
        if h in ('n_outliers',) and 'n_outliers' not in idx:
            idx['n_outliers'] = i
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--col-strs-id', default=None)
    ap.add_argument('--col-gene', default=None)
    ap.add_argument('--col-n-outliers', default=None)
    args = ap.parse_args()

    with open(args.catalog) as fh:
        r = csv.reader(fh, delimiter='\t')
        header = next(r)
        cols = detect_cols(header)
        for k, v in (('strs_id', args.col_strs_id), ('gene_name', args.col_gene),
                     ('n_outliers', args.col_n_outliers)):
            if v is not None:
                cols[k] = header.index(v) if v in header else int(v)
        sys.stderr.write(f"Colunas: {cols}\n")
        missing = [k for k, v in cols.items() if v is None]
        if missing:
            sys.exit(f"ERRO: colunas nao encontradas: {missing}")
        ci, cg, cn = cols['strs_id'], cols['gene_name'], cols['n_outliers']

        gene_loci = {}   # gene -> set(STRs_ID)
        nrows = 0
        for row in r:
            nrows += 1
            try:
                n = int(row[cn])
            except (ValueError, IndexError):
                continue
            if n <= 0:
                continue
            gene = row[cg] if cg < len(row) else ''
            sid = row[ci] if ci < len(row) else ''
            if not gene:
                continue
            gene_loci.setdefault(gene, set()).add(sid)

    sys.stderr.write(f"Catalogo: {nrows} linhas. Genes com >=1 STR outlier: {len(gene_loci)}\n")

    with open(args.out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene_name', 'n_outlier_str_loci'])
        for gene, s in sorted(gene_loci.items(), key=lambda kv: -len(kv[1])):
            w.writerow([gene, len(s)])
    sys.stderr.write(f"Escrevi {len(gene_loci)} genes -> {args.out}\n")


if __name__ == '__main__':
    main()
