#!/usr/bin/env python3
# annotate_catalog.py
# Enriquece os arquivos de 6.4.1.2 com a ANOTACAO COMPLETA do catalogo de STR da
# coorte (STRs_analysis_dataset.tsv), trazendo TODAS as colunas do catalogo
# EXCETO as relativas ao DBSCAN global (n_outliers, outlier_samples,
# outlier_residuals, n_clusters, noise_ratio).
#
# Colunas trazidas (ex.): group, age, sex, allele1_est, allele2_est, depth,
# repeat_unit, gene_id, gene_name, region, chrom, start, end, sample_id,
# pop (ancestralidade), contribution, type.
#
# O join e ADITIVO: preserva as colunas ja existentes no arquivo de entrada e
# so acrescenta as do catalogo que faltam. O ID (STRs_ID/strs_id) e o sample_id
# sao resolvidos de forma case-insensitive.
#
# Uso:
#   python3 annotate_catalog.py --catalog $CATALOG --in <arquivo> \
#       --id-col strs_id [--sample-col sample_id] --out <arquivo_anotado>
import argparse, csv, sys


DBSCAN_COLS = {'n_outliers', 'outlier_samples', 'outlier_residuals',
                'n_clusters', 'noise_ratio'}

# Colunas especificas de amostra: nao fazem sentido no nivel do STR (modo str).
SAMPLE_SPECIFIC = {'group', 'age', 'sex', 'allele1_est', 'allele2_est',
                    'depth', 'sample_id', 'pop', 'contribution'}


def get_ci(row, name):
    """Le coluna por nome case-insensitive (STRs_ID vs strs_id)."""
    ln = name.lower()
    for k, v in row.items():
        if k.lower() == ln:
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--in', required=True, dest='inp')
    ap.add_argument('--id-col', default='strs_id')
    ap.add_argument('--sample-col', default='sample_id')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    # ---- le catalogo e monta mapas ----
    sample_map = {}   # (id, sample) -> {col: val}
    str_map = {}      # id -> {col: val (primeiro), 'pops': set()}
    with open(args.catalog) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            sid = get_ci(row, 'STRs_ID') or ''
            samp = get_ci(row, 'sample_id') or ''
            sample_map[(sid, samp)] = row
            if sid not in str_map:
                str_map[sid] = {'_row': row, 'pops': set()}
            p = get_ci(row, 'pop') or ''
            if p:
                str_map[sid]['pops'].add(p)

    use_sample = False
    with open(args.inp) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        in_cols = list(r.fieldnames)
        # decide modo pelo header: se tem sample_id, join por (id, sample)
        if any(c.lower() == args.sample_col.lower() for c in in_cols):
            use_sample = True
        # colunas do catalogo a copiar (exceto dbscan e ja existentes)
        copy_cols = [c for c in r.fieldnames if c.lower() not in DBSCAN_COLS]
        # (r.fieldnames aqui sao do catalogo; precisamos da lista de colunas do catalogo)
    # reler catalogo p/ lista de colunas
    with open(args.catalog) as fh:
        r0 = csv.DictReader(fh, delimiter='\t')
        cat_cols = list(r0.fieldnames)
    copy_cols = [c for c in cat_cols if c.lower() not in DBSCAN_COLS]
    if not use_sample:
        # no nivel do STR: ignora colunas especificas de amostra
        copy_cols = [c for c in copy_cols if c.lower() not in SAMPLE_SPECIFIC]
    add_cols = [c for c in copy_cols if c not in in_cols]

    sys.stderr.write(f"Catalogo: {len(sample_map)} (id,sample); {len(str_map)} STRs unicos. "
                     f"Modo={'sample' if use_sample else 'str'}; colunas a acrescentar={len(add_cols)}\n")

    # Le tudo para memoria primeiro (seguro quando --in == --out, truncate antes de ler)
    with open(args.inp) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))

    n_ann = 0
    with open(args.out, 'w', newline='') as out:
        w = csv.DictWriter(out, delimiter='\t', fieldnames=in_cols + add_cols,
                           restval='', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            sid = get_ci(row, args.id_col) or ''
            cat = None
            if use_sample:
                samp = get_ci(row, args.sample_col) or ''
                cat = sample_map.get((sid, samp))
            if cat is None:
                rec = str_map.get(sid)
                if rec is not None:
                    cat = dict(rec['_row'])
                    cat['pops'] = ';'.join(sorted(rec['pops']))
            if cat:
                for c in add_cols:
                    row[c] = cat.get(c, '')
                n_ann += 1
            w.writerow(row)
    sys.stderr.write(f"Anotadas {n_ann} linhas -> {args.out}\n")


if __name__ == '__main__':
    main()
