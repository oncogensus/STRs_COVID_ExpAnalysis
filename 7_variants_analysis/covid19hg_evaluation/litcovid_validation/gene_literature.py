#!/usr/bin/env python3
# gene_literature.py
# Para cada gene outlier (data/outlier_genes.txt, da Parte 1), recupera os artigos
# COVID da literatura do LitCovid que mencionam o gene, usando o arquivo
# gene-anotado do LitCovid (litcovid2pubtator.json.gz, baixado por download_litcovid.sh).
#
# Mapeamento símbolo -> Entrez via MyGene.info; fallback por texto do mention quando
# o símbolo não mapeia. O arquivo PubTator é percorrido em streaming (JSON-lines, gzip),
# sem carregar tudo na memória.
#
# Saída:
#   results/gene_literature.tsv        (gene, entrez, pmid, title, journal, year)
#   results/gene_literature_summary.tsv (gene, n_outlier_str_loci, n_articles)
#
# Uso:
#   python3 gene_literature.py \
#       --genes data/outlier_genes.txt \
#       --litcovid data/litcovid2pubtator.json.gz \
#       --out results/gene_literature.tsv
import argparse, csv, gzip, json, sys, time, urllib.parse, urllib.request

MYGENE_URL = "https://mygene.info/v3/query"


def chunked(it, n):
    out = []
    for x in it:
        out.append(x)
        if len(out) == n:
            yield out
            out = []
    if out:
        yield out


def map_symbols_to_entrez(symbols):
    entrez_to_symbol = {}
    unmapped = set(symbols)
    for syms in chunked(symbols, 50):
        q = " OR ".join(f"symbol:{s}" for s in syms)
        params = urllib.parse.urlencode({'q': q, 'species': 'human',
                                         'fields': 'entrezgene,symbol'})
        url = f"{MYGENE_URL}?{params}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
        except Exception as e:
            sys.stderr.write(f"  mygene falhou para chunk: {e}\n")
            continue
        for hit in data.get('hits', []):
            ent = hit.get('entrezgene') or hit.get('_id')
            if ent is None:
                continue
            ent = str(ent)
            sym = hit.get('symbol') or ''
            entrez_to_symbol[ent] = sym
            if sym in unmapped:
                unmapped.discard(sym)
        time.sleep(0.1)
    return entrez_to_symbol, unmapped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genes', required=True)
    ap.add_argument('--litcovid', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--summary-out', default='results/gene_literature_summary.tsv')
    args = ap.parse_args()

    target_symbols = []
    gene_loci = {}
    with open(args.genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            g = row['gene_name']
            try:
                nl = int(row['n_outlier_str_loci'])
            except ValueError:
                nl = 0
            gene_loci[g] = nl
            target_symbols.append(g)
    sys.stderr.write(f"Genes-alvo: {len(target_symbols)}\n")

    sys.stderr.write("Mapeando simbolos -> Entrez (MyGene.info)...\n")
    entrez_to_symbol, unmapped = map_symbols_to_entrez(target_symbols)
    target_entrez = set(entrez_to_symbol.keys())
    symbol_to_entrez = {v: k for k, v in entrez_to_symbol.items()}
    sys.stderr.write(f"  mapeados: {len(target_entrez)}; "
                     f"nao mapeados (fallback texto): {len(unmapped)}\n")

    target_symbol_lower = {s.lower() for s in target_symbols}
    symbol_lower_to_symbol = {s.lower(): s for s in target_symbols}

    seen = set()
    articles = []
    n_total = 0
    with gzip.open(args.litcovid, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                art = json.loads(line)
            except Exception:
                continue
            n_total += 1
            if n_total % 500000 == 0:
                sys.stderr.write(f"  processados {n_total} artigos...\n")
            pmid = str(art.get('pmid', ''))
            title = art.get('title', '') or ''
            journal = art.get('journal', '') or ''
            year = art.get('publish_year', '') or ''
            matched = set()
            for ann in art.get('annotations', []):
                if ann.get('type') != 'Gene':
                    continue
                for m in ann.get('mentions', []):
                    gid = str(m.get('gene_id', ''))
                    if gid in target_entrez:
                        matched.add(entrez_to_symbol[gid])
                    else:
                        txt = (m.get('text', '') or '').lower()
                        if txt in target_symbol_lower:
                            matched.add(symbol_lower_to_symbol[txt])
            for g in matched:
                key = (g, pmid)
                if key in seen:
                    continue
                seen.add(key)
                articles.append((g, symbol_to_entrez.get(g, ''), pmid,
                                 title, journal, year))

    sys.stderr.write(f"Artigos LitCovid total: {n_total}; "
                     f"artigos associados aos genes-alvo: {len(articles)}\n")

    articles.sort(key=lambda x: (x[0], x[2]))
    with open(args.out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene', 'entrez', 'pmid', 'title', 'journal', 'year'])
        w.writerows(articles)

    counts = {}
    for g, _, _, _, _, _ in articles:
        counts[g] = counts.get(g, 0) + 1
    with open(args.summary_out, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['gene_name', 'n_outlier_str_loci', 'n_articles'])
        for g in target_symbols:
            w.writerow([g, gene_loci.get(g, 0), counts.get(g, 0)])

    sys.stderr.write(f"Escrevi {len(articles)} linhas -> {args.out}\n")
    sys.stderr.write(f"Escrevi resumo -> {args.summary_out}\n")


if __name__ == '__main__':
    main()
