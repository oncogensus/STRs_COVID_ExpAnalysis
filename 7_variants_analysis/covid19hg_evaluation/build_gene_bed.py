#!/usr/bin/env python3
# build_gene_bed.py
# Converte o GTF GENCODE v39 (hg38) em um BED de corpos gênicos:
#   chrom  start(0-based)  end(1-based)  gene_name
# Agrupa por gene_id pegando min(start)/max(end) (corpo do gene).
# Cromossomos normalizados para o formato "chrN"/"chrX".
#
# Uso:
#   python3 build_gene_bed.py \
#       --gtf data/gencode.v39.primary_assembly.annotation.gtf.gz \
#       --out genes.hg38.bed
import argparse, gzip, re, sys

def norm_chrom(raw):
    s = str(raw).strip()
    s = re.sub(r'^chr', '', s, flags=re.I)
    m = {'23': 'X', '24': 'Y', '25': 'MT', 'X': 'X', 'Y': 'Y', 'MT': 'MT', 'M': 'MT'}
    if s in m:
        s = m[s]
    return 'chr' + s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    gene_re = re.compile(r'gene_id "([^"]+)"')
    name_re = re.compile(r'gene_name "([^"]+)"')
    bodies = {}  # gene_id -> [chrom, min_start, max_end, gene_name]
    with gzip.open(args.gtf, 'rt') if args.gtf.endswith('.gz') else open(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene':
                continue
            chrom = norm_chrom(f[0])
            try:
                start = int(f[3]); end = int(f[4])
            except ValueError:
                continue
            gid = gene_re.search(f[8])
            nm = name_re.search(f[8])
            if not gid:
                continue
            gid = gid.group(1)
            name = nm.group(1) if nm else gid
            if gid not in bodies:
                bodies[gid] = [chrom, start, end, name]
            else:
                b = bodies[gid]
                if start < b[1]: b[1] = start
                if end > b[2]: b[2] = end

    rows = []
    for gid, b in bodies.items():
        rows.append((b[0], b[1], b[2], b[3]))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, 'w') as out:
        for chrom, start, end, name in rows:
            out.write(f"{chrom}\t{start-1}\t{end}\t{name}\n")  # BED 0-based start
    sys.stderr.write(f"Escrevi {len(rows)} genes -> {args.out}\n")

if __name__ == '__main__':
    main()
