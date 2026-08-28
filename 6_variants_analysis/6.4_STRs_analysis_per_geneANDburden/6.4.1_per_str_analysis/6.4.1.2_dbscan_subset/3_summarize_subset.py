#!/usr/bin/env python3
# summarize_subset.py
# Junta a tabela gene x STR sugestivo (overlap_suggestive_strs.py) com o resultado do
# DBSCAN re-rodado no subset (run_dbscan_subset.R), filtrando apenas os STRs que
# tiveram >=1 outlier (n_outliers>0).
#
# Uso:
#   python3 summarize_subset.py \
#       --overlap results/suggestive_gene_strs.tsv \
#       --dbscan results/suggestive_strs_outliers.tsv \
#       --out results/covid_suggestive_genes_with_outlier_STRs.tsv
import argparse, csv, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--overlap', required=True)
    ap.add_argument('--dbscan', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    dbscan = {}
    with open(args.dbscan) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            try:
                no = int(row['n_outliers'])
            except ValueError:
                no = 0
            dbscan[row['STRs_ID']] = (no, row.get('outlier_samples', ''),
                                      row.get('outlier_residuals', ''),
                                      row.get('n_clusters', ''),
                                      row.get('noise_ratio', ''),
                                      row.get('eps', ''),
                                      row.get('minPts', ''),
                                      row.get('cutoff', ''),
                                      row.get('max_residual', ''),
                                      row.get('mean_residual', ''))

    n_out = 0
    with open(args.overlap) as fh, open(args.out, 'w', newline='') as out:
        r = csv.DictReader(fh, delimiter='\t')
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene', 'chrom', 'gene_start', 'gene_end', 'best_p', 'phenotypes',
                    'lead_snp', 'strs_id', 'str_chrom', 'str_start', 'str_end',
                    'repeatunit', 'n_carriers', 'n_outliers', 'outlier_samples',
                    'outlier_residuals', 'n_clusters', 'noise_ratio', 'eps', 'minPts',
                    'cutoff', 'max_residual', 'mean_residual'])
        for row in r:
            sid = row['strs_id']
            if sid not in dbscan:
                continue
            no, os_, ore, nc, nr, eps, mpts, cut, mres, meres = dbscan[sid]
            if no <= 0:
                continue
            w.writerow([row['gene'], row['chrom'], row['gene_start'], row['gene_end'],
                        row['best_p'], row['phenotypes'], row['lead_snp'], sid,
                        row['str_chrom'], row['str_start'], row['str_end'],
                        row['repeatunit'], row['n_carriers'], no, os_, ore, nc, nr,
                        eps, mpts, cut, mres, meres])
            n_out += 1
    sys.stderr.write(f"Escrevi {n_out} pares (gene sugestivo x STR outlier) -> {args.out}\n")


if __name__ == '__main__':
    main()
