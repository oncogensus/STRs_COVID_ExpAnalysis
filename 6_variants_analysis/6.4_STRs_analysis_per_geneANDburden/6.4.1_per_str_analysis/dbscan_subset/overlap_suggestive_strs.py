#!/usr/bin/env python3
# overlap_suggestive_strs.py
# Cruzamento coord-a-coord:
#   genes COVID sugestivos (covid_genes_suggestive.tsv, p<1e-5) x catálogo de STRs
#   da coorte. Para cada gene sugestivo, lista os loci de STR da coorte que caem no
#   corpo do gene +/- window. Gera a lista de STRs_ID que serão submetidos ao DBSCAN
#   (run_dbscan_subset.R) e a tabela gene x STR.
#
# Uso:
#   python3 overlap_suggestive_strs.py \
#       --catalog $CATALOG \
#       --covid-genes results/covid_genes_suggestive.tsv \
#       --out results/suggestive_gene_strs.tsv \
#       --ids-out data/suggestive_strs_ids.txt
import argparse, csv, re, sys


def norm_chrom(raw):
    s = str(raw).strip()
    s = re.sub(r'^chr', '', s, flags=re.I)
    m = {'23': 'X', '24': 'Y', '25': 'MT', 'X': 'X', 'Y': 'Y', 'MT': 'MT', 'M': 'MT'}
    if s in m:
        s = m[s]
    return 'chr' + s


def detect_cols(header):
    low = [h.lower().lstrip('#') for h in header]
    idx = {k: None for k in ('chrom', 'left', 'right', 'motif', 'sample', 'strs_id')}
    for i, h in enumerate(low):
        if idx['chrom'] is None and h == 'chrom':
            idx['chrom'] = i
        if idx['left'] is None and ('left' in h or h == 'start' or 'start' in h or 'begin' in h):
            idx['left'] = i
        if idx['right'] is None and ('right' in h or h == 'end' or 'end' in h):
            idx['right'] = i
        if idx['motif'] is None and ('repeat' in h or 'motif' in h):
            idx['motif'] = i
        if idx['sample'] is None and 'sample' in h:
            idx['sample'] = i
        if idx['strs_id'] is None and 'strs_id' in h:
            idx['strs_id'] = i
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--covid-genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ids-out', required=True)
    ap.add_argument('--window', type=int, default=50000)
    ap.add_argument('--col-chrom', default=None)
    ap.add_argument('--col-left', default=None)
    ap.add_argument('--col-right', default=None)
    ap.add_argument('--col-motif', default=None)
    ap.add_argument('--col-sample', default=None)
    ap.add_argument('--col-strs-id', default=None)
    args = ap.parse_args()

    genes = {}
    with open(args.covid_genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            chrom = norm_chrom(row['chrom'])
            try:
                gs = int(row['gene_start'])
                ge = int(row['gene_end'])
            except ValueError:
                continue
            genes.setdefault(chrom, []).append(
                [gs - args.window, ge + args.window, row['gene'],
                 row['best_p'], row['phenotypes'], row['lead_snp']])

    with open(args.catalog) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        cols = detect_cols(header)
        for key, val in (('chrom', args.col_chrom), ('left', args.col_left),
                         ('right', args.col_right), ('motif', args.col_motif),
                         ('sample', args.col_sample), ('strs_id', args.col_strs_id)):
            if val is not None:
                cols[key] = header.index(val) if val in header else int(val)
        sys.stderr.write(f"Colunas detectadas: {cols}\n")
        missing = [k for k, v in cols.items() if v is None]
        if missing:
            sys.exit(f"ERRO: colunas nao encontradas: {missing}")
        ci = cols['strs_id']

        loci = {}
        nrows = 0
        for row in reader:
            nrows += 1
            try:
                chrom = norm_chrom(row[cols['chrom']])
                l = int(row[cols['left']])
                rg = int(row[cols['right']])
            except (ValueError, IndexError):
                continue
            motif = row[cols['motif']] if cols['motif'] < len(row) else ''
            sample = row[cols['sample']] if cols['sample'] < len(row) else ''
            s = min(l, rg)
            e = max(l, rg)
            sid = row[ci] if (ci is not None and ci < len(row)) else f"{chrom}:{l}-{rg}:{motif}"
            key = (chrom, l, rg, motif)
            if key not in loci:
                loci[key] = [chrom, s, e, motif, sid, set()]
            loci[key][5].add(sample)

    sys.stderr.write(f"Catalogo: {nrows} linhas, {len(loci)} loci.\n")

    out_rows = []
    idset = set()
    for key, (chrom, s, e, motif, sid, samples) in loci.items():
        for g in genes.get(chrom, []):
            if g[0] <= e and g[1] >= s:
                out_rows.append([g[2], chrom, g[0] + args.window, g[1] - args.window,
                                 g[3], g[4], g[5], sid, chrom, s, e, motif, len(samples)])
                idset.add(sid)
    out_rows.sort(key=lambda x: (x[0], float(x[4]) if x[4] != '' else 1.0))

    with open(args.out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene', 'chrom', 'gene_start', 'gene_end', 'best_p', 'phenotypes',
                    'lead_snp', 'strs_id', 'str_chrom', 'str_start', 'str_end',
                    'repeatunit', 'n_carriers'])
        w.writerows(out_rows)
    with open(args.ids_out, 'w') as fid:
        for sid in sorted(idset):
            fid.write(sid + "\n")

    genes_hit = sorted({r[0] for r in out_rows})
    sys.stderr.write(f"Escrevi {len(out_rows)} pares -> {args.out}; "
                     f"{len(idset)} STRs_ID unicos -> {args.ids_out}\n")
    sys.stderr.write(f"Genes sugestivos com STR na coorte: {len(genes_hit)}\n")
    if genes_hit:
        sys.stderr.write("  " + ", ".join(genes_hit) + "\n")


if __name__ == '__main__':
    main()
