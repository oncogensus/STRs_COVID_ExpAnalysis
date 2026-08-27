#!/usr/bin/env python3
# extract_covid_genes.py
# Le os summary stats do COVID-19 HG r7 (A2/B2/C2, leave_23andme), extrai os SNPs
# com p < --p-thresh (default 5e-8) e mapeia cada um ao gene se cai no corpo
# do gene ou ate --window bp flanqueando (default 50000). Gera um conjunto de
# "genes COVID" com o melhor p por gene e o fenotipo(s).
#
# Uso:
#   python3 extract_covid_genes.py \
#       --sumstats data/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz \
#                  data/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz \
#                  data/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz \
#       --gene-bed genes.hg38.bed \
#       --out covid_genes.tsv
import argparse, gzip, re, sys, os, bisect

def norm_chrom(raw):
    s = str(raw).strip()
    s = re.sub(r'^chr', '', s, flags=re.I)
    m = {'23': 'X', '24': 'Y', '25': 'MT', 'X': 'X', 'Y': 'Y', 'MT': 'MT', 'M': 'MT'}
    if s in m:
        s = m[s]
    return 'chr' + s

def phenotype_of(path):
    base = os.path.basename(path)
    m = re.search(r'_(A\d|B\d|C\d)_', base)
    return m.group(1) if m else base

def load_genes(path, window):
    by_chr = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            chrom = f[0]
            try:
                s = int(f[1]); e = int(f[2])
            except ValueError:
                continue
            name = f[3]
            by_chr.setdefault(chrom, []).append([s, e, name])
    for c in by_chr:
        by_chr[c].sort(key=lambda x: x[0])
    return by_chr

def map_snp(by_chr, chrom, pos, window):
    genes = by_chr.get(chrom)
    if not genes:
        return []
    # candidatos: start <= pos+window (depois filtra por end+window >= pos)
    idx = bisect.bisect_right([g[0] for g in genes], pos + window)
    hits = []
    for g in genes[:idx]:
        if g[1] + window >= pos:
            hits.append((g[2], g[0], g[1]))  # (name, gstart, gend)
    return hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sumstats', nargs='+', required=True)
    ap.add_argument('--gene-bed', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--window', type=int, default=50000)
    ap.add_argument('--p-thresh', type=float, default=5e-8)
    args = ap.parse_args()

    sys.stderr.write("Carregando genes...\n")
    by_chr = load_genes(args.gene_bed, args.window)

    # gene -> [best_p, set(phenotypes), lead_snp, chrom, gstart, gend]
    genes = {}

    for path in args.sumstats:
        pheno = phenotype_of(path)
        sys.stderr.write(f"Processando {path} (fenotipo {pheno})...\n")
        fh = gzip.open(path, 'rt') if path.endswith('.gz') else open(path)
        header = fh.readline().rstrip('\n').lstrip('#').split('\t')
        col = {h: i for i, h in enumerate(header)}
        i_chr = col.get('CHR')
        i_pos = col.get('POS')
        i_p = col.get('all_inv_var_meta_p')
        i_beta = col.get('all_inv_var_meta_beta')
        i_rs = col.get('rsid')
        if i_rs is None:
            i_rs = col.get('SNP')
        if i_chr is None or i_pos is None or i_p is None:
            sys.exit(f"ERRO: colunas esperadas nao encontradas em {path}. Header: {header}")
        n = 0
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= max(i_chr, i_pos, i_p):
                continue
            try:
                p = float(f[i_p])
            except ValueError:
                continue
            if p <= 0 or p >= args.p_thresh:
                continue
            try:
                pos = int(f[i_pos])
            except ValueError:
                continue
            chrom = norm_chrom(f[i_chr])
            rs = f[i_rs] if (i_rs is not None and i_rs < len(f) and f[i_rs]) else f"{chrom}:{pos}"
            for gname, gs, ge in map_snp(by_chr, chrom, pos, args.window):
                if gname not in genes:
                    genes[gname] = [p, {pheno}, rs, chrom, gs, ge]
                else:
                    rec = genes[gname]
                    rec[1].add(pheno)
                    if p < rec[0]:
                        rec[0] = p
                        rec[2] = rs
            n += 1
        fh.close()
        sys.stderr.write(f"  {n} SNPs significativos (p<{args.p_thresh}) em {pheno}\n")

    with open(args.out, 'w') as out:
        out.write("gene\tchrom\tgene_start\tgene_end\tbest_p\tphenotypes\tlead_snp\n")
        for gname, rec in sorted(genes.items(), key=lambda kv: kv[1][0]):
            # ajusta gene_start/end para o corpo do gene (do primeiro SNP nao eh ideal;
            # preenche abaixo com o corpo real quando disponivel)
            out.write(f"{gname}\t{rec[3]}\t{rec[4]}\t{rec[5]}\t{rec[0]:.3e}\t{','.join(sorted(rec[1]))}\t{rec[2]}\n")
    sys.stderr.write(f"Escrevi {len(genes)} genes COVID -> {args.out}\n")

if __name__ == '__main__':
    main()
