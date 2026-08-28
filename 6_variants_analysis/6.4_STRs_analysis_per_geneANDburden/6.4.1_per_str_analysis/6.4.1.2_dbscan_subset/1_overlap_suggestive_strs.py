#!/usr/bin/env python3
# overlap_suggestive_strs.py
# ---------------------------------------------------------------------------
# PROPOSITO
#   Cruzamento coord-a-coord (genoma hg38) entre:
#     (A) genes COVID-19 sugestivos  -> covid_genes_suggestive.tsv (p < 1e-5,
#         fenotipos A2/B2/C2 do COVID-19 HG r7, leave_23andme), e
#     (B) catalogo completo de STRs da coorte (STRling) -> $CATALOG.
#
#   Para cada gene sugestivo, lista os LOCI de STR da coorte que caem NO CORPO
#   DO GENE. O resultado alimenta o DBSCAN subset
#   (2_run_dbscan_subset.R / .sh): so os STRs_ID desta lista sao re-analisados.
#
# ALGORITMO (resumo)
#   1. Carrega genes sugestivos -> dicionario {chrom: [(gs, ge, ...)]}
#   2. Carrega catalogo de STR -> agrega por locus unico (chrom,left,right,motif),
#      acumulando o conjunto de amostras (carriers) de cada locus.
#   3. Intersecao de intervalos: gene_overlap se (gene_start <= str_end) e
#      (gene_end >= str_start). Coleta (gene, strs_id) para cada par sobreposto
#      (STR que toca o corpo do gene, gene_start..gene_end).
#   4. Escreve tabela gene x STR e o arquivo de IDs (um STRs_ID por linha).
#
# ENTRADAS
#   --catalog        Catalogo de STR da coorte (STRs_analysis_dataset.tsv):
#                    colunas cromossomo, left, right, repeatunit, sample, strs_id.
#   --covid-genes    Saida de extract_covid_genes.py filtrada a p<1e-5
#                    (results/covid_genes_suggestive.tsv): chrom, gene_start,
#                    gene_end, gene, best_p, phenotypes, lead_snp.
#   --out            Tabela de saida gene x STR (results/suggestive_gene_strs.tsv).
#   --ids-out        Lista de STRs_ID unicos (data/suggestive_strs_ids.txt) usada
#                    pelo subset do DBSCAN e pelo crossvalidate LitCovid.
#   --debug          Ativa mensagens de debug verbosas (progresso, contagens).
#
# SAIDAS
#   --out    12 colunas: gene, chrom, gene_start, gene_end, best_p, phenotypes,
#            lead_snp, strs_id, str_chrom, str_start, str_end, repeatunit, n_carriers
#   --ids-out  Um STRs_ID por linha (deduplicado).
# ---------------------------------------------------------------------------
# Uso:
#   python3 overlap_suggestive_strs.py \
#       --catalog $CATALOG \
#       --covid-genes results/covid_genes_suggestive.tsv \
#       --out results/suggestive_gene_strs.tsv \
#       --ids-out data/suggestive_strs_ids.txt
# ---------------------------------------------------------------------------
import argparse, csv, re, sys, time


def dbg(msg, debug=False):
    """Escreve msg em stderr apenas se debug=True. Padrao: silencioso."""
    if debug:
        sys.stderr.write(f"[DEBUG] {msg}\n")


def norm_chrom(raw):
    """Normaliza qualquer rotulo de cromossomo para o formato 'chrN'/'chrX'...

    Aceita tanto 'chr3' quanto '3', e mapeia aliases numericos usados por alguns
    summary stats (23->X, 24->Y, 25/MT/M->MT). Garante que genes e STRs comparados
    usem a MESMA nomenclatura de cromossomo antes do cruzamento coord-a-coord.
    """
    s = str(raw).strip()
    s = re.sub(r'^chr', '', s, flags=re.I)          # remove prefixo 'chr' se houver
    m = {'23': 'X', '24': 'Y', '25': 'MT', 'X': 'X', 'Y': 'Y', 'MT': 'MT', 'M': 'MT'}
    if s in m:
        s = m[s]
    return 'chr' + s


def detect_cols(header):
    """Infere os indices das colunas essenciais no catalogo de STR por heuristica.

    Procura nomes canonicos (chrom, left/start, right/end, repeat/motif, sample,
    strs_id). Retorna dicionario {nome: indice|None}. Indices None indicam coluna
    nao encontrada automaticamente (podem ser sobrescritos via --col-*).
    """
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
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--covid-genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ids-out', required=True)
    ap.add_argument('--col-chrom', default=None)
    ap.add_argument('--col-left', default=None)
    ap.add_argument('--col-right', default=None)
    ap.add_argument('--col-motif', default=None)
    ap.add_argument('--col-sample', default=None)
    ap.add_argument('--col-strs-id', default=None)
    ap.add_argument('--debug', action='store_true',
                    help='Ativa mensagens de debug verbosas em stderr')
    args = ap.parse_args()
    D = args.debug

    dbg(f"args: catalog={args.catalog} covid-genes={args.covid_genes} "
        f"out={args.out} ids-out={args.ids_out}", D)

    # ----------------------------------------------------------------------
    # 1) Carregar genes COVID sugestivos
    #    genes: {chrom_normalizado: [(gs, ge, gene, best_p,
    #                                 phenotypes, lead_snp), ...]}
    #    Usamos as coordenadas EXATAS do corpo do gene (gene_start, gene_end).
    #    So STR que sobrepoe esse intervalo sera considerado.
    # ----------------------------------------------------------------------
    genes = {}
    with open(args.covid_genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        n_genes_raw = 0
        for row in r:
            n_genes_raw += 1
            chrom = norm_chrom(row['chrom'])
            try:
                gs = int(row['gene_start'])
                ge = int(row['gene_end'])
            except ValueError:
                # linha com coordenada invalida -> pula silenciosamente
                continue
            genes.setdefault(chrom, []).append(
                [gs, ge, row['gene'],
                 row['best_p'], row['phenotypes'], row['lead_snp']])
    n_genes = sum(len(v) for v in genes.values())
    sys.stderr.write(f"Genes sugestivos: {n_genes} validos de {n_genes_raw} linhas "
                     f"(em {len(genes)} cromossomos)\n")
    dbg("Genes por cromossomo: " + ", ".join(
        f"{c}={len(v)}" for c, v in sorted(genes.items())), D)

    # ----------------------------------------------------------------------
    # 2) Carregar catalogo de STR e agregar por locus unico
    #    loci: {(chrom, left, right, motif): [chrom, s, e, motif, strs_id,
    #                                         set(samples)]}
    #    - Chave de 4 campos dedupica loci identicos (mesmo repeat, mesma posicao)
    #      independente da amostra.
    #    - s/e = min/max de left/right (garante inicio<=fim mesmo se trocados).
    #    - O set de samples conta quantas amostras-carrier tem aquele locus.
    # ----------------------------------------------------------------------
    with open(args.catalog) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        cols = detect_cols(header)
        # sobrescreve colunas detectadas automaticamente com os --col-* explicitos
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
        n_skip = 0
        for row in reader:
            nrows += 1
            try:
                chrom = norm_chrom(row[cols['chrom']])
                l = int(row[cols['left']])
                rg = int(row[cols['right']])
            except (ValueError, IndexError):
                # coordenada ausente/ilegivel -> pula a linha
                n_skip += 1
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
            # progresso de debug a cada 200k linhas do catalogo
            if D and nrows % 200000 == 0:
                dbg(f"catalogo: {nrows} linhas lidas, {len(loci)} loci ate agora", D)

    sys.stderr.write(f"Catalogo: {nrows} linhas ({n_skip} puladas), "
                     f"{len(loci)} loci unicos.\n")

    # ----------------------------------------------------------------------
    # 3) Intersecao de intervalos gene x STR
    #    Condicao de sobreposicao: gene_start <= str_end  E  gene_end >= str_start
    #    (intervalos fechados). g[0]/g[1] sao as coordenadas EXATAS do corpo do
    #    gene, entao so conta STR que toca o gene (gene_start..gene_end).
    # ----------------------------------------------------------------------
    out_rows = []
    idset = set()
    n_matches = 0
    for key, (chrom, s, e, motif, sid, samples) in loci.items():
        for g in genes.get(chrom, []):
            if g[0] <= e and g[1] >= s:               # sobreposicao gene<->STR
                out_rows.append([g[2], chrom, g[0], g[1],
                                 g[3], g[4], g[5], sid, chrom, s, e, motif, len(samples)])
                idset.add(sid)
                n_matches += 1
    dbg(f"Pares (gene,STR) sobrepostos encontrados: {n_matches}", D)
    out_rows.sort(key=lambda x: (x[0], float(x[4]) if x[4] != '' else 1.0))

    # ----------------------------------------------------------------------
    # 4) Escrever saidas (tabela gene x STR e lista de IDs)
    # ----------------------------------------------------------------------
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
    sys.stderr.write(f"Tempo total: {time.time() - t0:.1f}s\n")


if __name__ == '__main__':
    main()
