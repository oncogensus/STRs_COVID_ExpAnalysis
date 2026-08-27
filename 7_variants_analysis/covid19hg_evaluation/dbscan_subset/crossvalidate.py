#!/usr/bin/env python3
# crossvalidate.py
# Cruza os dois lados da validacao:
#   - LitCovid (Parte 1): genes outlier do DBSCAN global que tem literatura COVID
#     (results/gene_literature_summary.tsv, em litcovid_validation)
#   - DBSCAN subset (Parte 2): genes sugestivos do COVID-19 HG que contem STR outlier
#     na coorte (results/covid_suggestive_genes_with_outlier_STRs.tsv)
# Saida: results/crossvalidated_genes.tsv (genes em AMBOS os conjuntos) + resumo.
#
# Uso:
#   python3 crossvalidate.py \
#       --litcovid-summary ../litcovid_validation/results/gene_literature_summary.tsv \
#       --subset results/covid_suggestive_genes_with_outlier_STRs.tsv \
#       --out results/crossvalidated_genes.tsv
import argparse, csv, sys


ORIGINAL_8 = ["ROBO2", "ANK3", "CDH12", "NKAIN2", "SEMA6D",
              "KCNH1", "KCNQ5", "ST6GALNAC3"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--litcovid-summary', required=True)
    ap.add_argument('--subset', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--original8', default=",".join(ORIGINAL_8))
    args = ap.parse_args()
    original8 = set(g for g in args.original8.split(',') if g)

    lit = {}
    with open(args.litcovid_summary) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            lit[row['gene_name']] = (row.get('n_outlier_str_loci', ''),
                                     row.get('n_articles', ''))

    sub = {}
    with open(args.subset) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            g = row['gene']
            rec = sub.setdefault(g, {'best_p': row['best_p'],
                                     'phenotypes': row['phenotypes'],
                                     'strs': set(), 'n_outliers': 0})
            rec['strs'].add(row['strs_id'])
            try:
                rec['n_outliers'] += int(row['n_outliers'])
            except ValueError:
                pass

    inter = sorted(set(lit) & set(sub))
    sys.stderr.write(f"LitCovid genes: {len(lit)}; Subset genes: {len(sub)}; "
                     f"Intersecao: {len(inter)}\n")

    with open(args.out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene', 'n_outlier_str_loci', 'n_litcovid_articles',
                    'best_p_suggestive', 'phenotypes', 'n_outlier_strs_in_gene',
                    'n_outliers_total', 'in_original_8'])
        for g in inter:
            lo, na = lit[g]
            rec = sub[g]
            w.writerow([g, lo, na, rec['best_p'], rec['phenotypes'],
                        len(rec['strs']), rec['n_outliers'],
                        'yes' if g in original8 else 'no'])
    sys.stderr.write(f"Escrevi {len(inter)} genes cross-validados -> {args.out}\n")
    if inter:
        sys.stderr.write("  " + ", ".join(inter) + "\n")


if __name__ == '__main__':
    main()
