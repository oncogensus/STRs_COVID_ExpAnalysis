#!/usr/bin/env python3
# annotate_location.py
# Anexa a LOCALIZACAO GENOMICA do STR (str_chrom, str_start, str_end, repeatunit,
# str_locus) a um TSV que contenha a coluna STRs_ID, fazendo join com a tabela
# de overlap (suggestive_gene_strs.tsv), que ja traz as coordenadas do locus.
#
# Usado em 6.4.1.2 para enriquecer os arquivos de residuos e de outliers do
# DBSCAN com a anotacao de onde cada STR esta no genoma (hg38).
#
# Uso:
#   python3 annotate_location.py \
#       --in  data/suggestive_strs_residuals.tsv \
#       --overlap results/suggestive_gene_strs.tsv \
#       --out data/suggestive_strs_residuals.tsv
#
# Se --out for igual a --in, escreve num temporario e substitui (in-place).
import argparse, csv, sys, os


def get_ci(row, name):
    """Le valor de coluna com nome case-insensitive (ex.: STRs_ID vs strs_id)."""
    ln = name.lower()
    for k, v in row.items():
        if k.lower() == ln:
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', required=True, dest='inp')
    ap.add_argument('--overlap', required=True,
                    help='Tabela gene x STR (suggestive_gene_strs.tsv) de onde vem as coordenadas.')
    ap.add_argument('--id-col', default='strs_id')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    # mapa STRs_ID -> (str_chrom, str_start, str_end, repeatunit); dedup por 1a ocorrencia
    loc = {}
    with open(args.overlap) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            sid = get_ci(row, args.id_col)
            if not sid or sid in loc:
                continue
            loc[sid] = (get_ci(row, 'str_chrom') or '', get_ci(row, 'str_start') or '',
                        get_ci(row, 'str_end') or '', get_ci(row, 'repeatunit') or '')

    sys.stderr.write(f"Mapa de localizacao: {len(loc)} STRs_ID unicos (de {args.overlap})\n")

    tmp = args.out
    if tmp == args.inp:
        tmp = args.out + '.annot.tmp'

    n_ann = 0
    with open(args.inp) as fh, open(tmp, 'w', newline='') as out:
        r = csv.DictReader(fh, delimiter='\t')
        w = csv.DictWriter(out, delimiter='\t', fieldnames=list(r.fieldnames) +
                           ['str_chrom', 'str_start', 'str_end', 'repeatunit', 'str_locus'])
        w.writeheader()
        for row in r:
            sid = get_ci(row, args.id_col) or ''
            c = loc.get(sid, ('', '', '', ''))
            row['str_chrom'] = c[0]
            row['str_start'] = c[1]
            row['str_end'] = c[2]
            row['repeatunit'] = c[3]
            row['str_locus'] = (f"{c[0]}:{c[1]}-{c[2]}:{c[3]}"
                                 if c[0] and c[1] and c[2] else '')
            if c[0]:
                n_ann += 1
            w.writerow(row)

    if tmp != args.out:
        os.replace(tmp, args.out)
    sys.stderr.write(f"Anotadas {n_ann} linhas com localizacao genomica -> {args.out}\n")


if __name__ == '__main__':
    main()
