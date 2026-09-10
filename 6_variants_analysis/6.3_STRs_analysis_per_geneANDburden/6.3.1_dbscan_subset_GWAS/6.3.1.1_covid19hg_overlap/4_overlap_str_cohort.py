#!/usr/bin/env python3
# overlap_str_cohort.py
# Cruzamento coord-a-coord:
#   genes COVID (covid_genes.tsv)  x  catalogo completo de STRs da coorte
#   (STRs_analysis_dataset.tsv / global_STRs_filtered.tsv).
# Um STR "cai no gene" se sua extensao [s,e] sobrepoe o corpo do gene
# (gene_start..gene_end), SEM janela de flanco.
#
# Uso:
#   python3 overlap_str_cohort.py \
#       --catalog /storage2/.../samples/STRs_analysis_dataset.tsv \
#       --covid-genes covid_genes.tsv \
#       --out covid_genes_with_cohort_STRs.tsv
#       [--debug]
#
# DEBUG: com --debug, imprime no stderr o progresso do cruzamento, o numero
#   de amostras identificadas (unicas e totais) e o detalhamento por
#   localizacao (cromossomo) das STRs que caem em genes COVID.
#
# Colunas do catalogo sao auto-detectadas (chrom/left/right/repeatunit/sample);
# sobrescreva com --col-chrom/--col-left/--col-right/--col-motif/--col-sample se preciso.
import argparse, re, sys, csv

def norm_chrom(raw):
    s = str(raw).strip()
    s = re.sub(r'^chr', '', s, flags=re.I)
    m = {'23': 'X', '24': 'Y', '25': 'MT', 'X': 'X', 'Y': 'Y', 'MT': 'MT', 'M': 'MT'}
    if s in m:
        s = m[s]
    return 'chr' + s

def detect_cols(header):
    low = [h.lower().lstrip('#') for h in header]
    idx = {k: None for k in ('chrom', 'left', 'right', 'motif', 'sample')}
    for i, h in enumerate(low):
        if idx['chrom'] is None and h in ('chrom',):
            idx['chrom'] = i
        if idx['left'] is None and ('left' in h or h == 'start' or 'start' in h or 'begin' in h):
            idx['left'] = i
        if idx['right'] is None and ('right' in h or h == 'end' or 'end' in h):
            idx['right'] = i
        if idx['motif'] is None and ('repeat' in h or 'motif' in h):
            idx['motif'] = i
        if idx['sample'] is None and 'sample' in h:
            idx['sample'] = i
    return idx

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--covid-genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--col-chrom', default=None)
    ap.add_argument('--col-left', default=None)
    ap.add_argument('--col-right', default=None)
    ap.add_argument('--col-motif', default=None)
    ap.add_argument('--col-sample', default=None)
    ap.add_argument('--debug', action='store_true',
                   help='Imprime debugs de progresso, amostras identificadas '
                        'e detalhamento por localizacao (cromossomo).')
    args = ap.parse_args()

    def dbg(msg):
        if args.debug:
            sys.stderr.write(f"[DEBUG] {msg}\n")

    # ---- carrega genes COVID (corpo do gene, sem janela) ----
    covid = {}  # chrom -> list of [gstart, gend, name, best_p, pheno, lead]
    with open(args.covid_genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            chrom = norm_chrom(row['chrom'])
            try:
                gs = int(row['gene_start']); ge = int(row['gene_end'])
            except ValueError:
                continue
            covid.setdefault(chrom, []).append(
                [gs, ge, row['gene'], row['best_p'], row['phenotypes'], row['lead_snp']])

    n_genes = sum(len(v) for v in covid.values())
    dbg(f"Genes COVID carregados: {n_genes} em {len(covid)} cromossomos.")
    if args.debug:
        per_chrom = ", ".join(f"{c}={len(v)}" for c, v in sorted(covid.items()))
        dbg(f"Genes COVID por cromossomo: {per_chrom}")
    sys.stderr.write(f"Genes COVID: {n_genes} (corpo do gene, sem janela)\n")

    # ---- le catalogo de STR ----
    with open(args.catalog) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        cols = detect_cols(header)
        # overrides
        for key, val in (('chrom', args.col_chrom), ('left', args.col_left),
                         ('right', args.col_right), ('motif', args.col_motif),
                         ('sample', args.col_sample)):
            if val is not None:
                cols[key] = header.index(val) if val in header else int(val)
        sys.stderr.write(f"Colunas detectadas: {cols}\n")
        missing = [k for k, v in cols.items() if v is None]
        if missing:
            sys.exit(f"ERRO: nao encontrei colunas {missing}. Use --col-* para especificar.")

        loci = {}  # key -> [chrom, s, e, motif, set(samples)]
        nrows = 0
        for row in reader:
            nrows += 1
            try:
                chrom = norm_chrom(row[cols['chrom']])
                l = int(row[cols['left']]); rg = int(row[cols['right']])
            except (ValueError, IndexError):
                continue
            motif = row[cols['motif']] if cols['motif'] < len(row) else ''
            sample = row[cols['sample']] if cols['sample'] < len(row) else ''
            s = min(l, rg); e = max(l, rg)
            key = (chrom, l, rg, motif)
            if key not in loci:
                loci[key] = [chrom, s, e, motif, set()]
            if sample:
                loci[key][4].add(sample)

    sys.stderr.write(f"Catalogo: {nrows} linhas, {len(loci)} loci unicos de STR.\n")
    dbg(f"Iniciando cruzamento coord-a-coord de {len(loci)} loci STR contra genes COVID...")

    # ---- cruzamento (somente dentro do corpo do gene) ----
    out_rows = []
    # acompanhamento para estatisticas de "amostras identificadas" e localizacao
    all_samples = set()          # amostras unicas que carregam algum STR no gene
    loci_hit = set()             # loci STR unicos que caem em algum gene COVID
    carr_by_chrom = {}           # chrom -> [n_pares, n_loci_str_unicos, n_amostras_unicas]
    for n_done, (key, (chrom, s, e, motif, samples)) in enumerate(loci.items(), 1):
        if args.debug and n_done % 100000 == 0:
            dbg(f"  processados {n_done}/{len(loci)} loci STR...")
        genes = covid.get(chrom, [])
        # candidatos: gstart <= e  (depois filtra gend >= s)
        for g in genes:
            if g[0] <= e and g[1] >= s:
                out_rows.append([
                    g[2], chrom, g[0], g[1],
                    g[3], g[4], g[5],
                    f"{chrom}:{s}-{e}:{motif}", chrom, s, e, motif, len(samples)
                ])
                loci_hit.add(key)
                all_samples |= samples
                st = carr_by_chrom.setdefault(chrom, [0, set(), set()])
                st[0] += 1
                st[1].add(key)
                st[2] |= samples

    out_rows.sort(key=lambda x: (x[0], float(x[4]) if x[4] != '' else 1.0))
    with open(args.out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene', 'chrom', 'gene_start', 'gene_end', 'best_p',
                    'phenotypes', 'lead_snp', 'str_locus', 'str_chrom',
                    'str_start', 'str_end', 'repeatunit', 'n_carriers'])
        w.writerows(out_rows)

    # ---- resumo / debug de amostras identificadas e localizacao ----
    total_carriers = sum(r[12] for r in out_rows)   # soma de ocorrencias (mesma amostra pode contar em >1 STR)
    genes_hit = sorted({r[0] for r in out_rows})
    sys.stderr.write(f"Escrevi {len(out_rows)} pares (gene COVID x STR coorte) -> {args.out}\n")
    sys.stderr.write(f"Genes COVID que TEM STR na coorte: {len(genes_hit)}\n")
    if genes_hit:
        sys.stderr.write("  " + ", ".join(genes_hit) + "\n")
    sys.stderr.write(
        f"AMOSTRAS IDENTIFICADAS: {len(all_samples)} unicas "
        f"(somando portadores={total_carriers}) em {len(loci_hit)} loci STR "
        f"distintos que caem em genes COVID.\n")

    if args.debug:
        dbg("Detalhamento por localizacao (cromossomo) das STRs que caem em genes COVID:")
        for chrom in sorted(carr_by_chrom):
            n_pairs, loci_set, samp_set = carr_by_chrom[chrom]
            dbg(f"  {chrom}: {n_pairs} pares, {len(loci_set)} loci STR, "
                f"{len(samp_set)} amostras unicas")
        # top genes por numero de loci STR distintos
        gene_to_loci = {}
        for r in out_rows:
            gene_to_loci.setdefault(r[0], set()).add((r[8], r[9], r[10], r[11]))
        dbg("Top genes por nº de loci STR distintos que caem no gene:")
        for gene, ls in sorted(gene_to_loci.items(), key=lambda kv: -len(kv[1]))[:10]:
            dbg(f"  {gene}: {len(ls)} loci STR")

if __name__ == '__main__':
    main()
