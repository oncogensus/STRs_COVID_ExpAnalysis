#!/usr/bin/env python3
# overlap_suggestive_strs.py
# ---------------------------------------------------------------------------
# PROPOSITO
#   Cruzamento por gene_name entre:
#     (A) genes COVID-19 sugestivos  -> covid_genes_suggestive.tsv (p < 1e-5,
#         fenotipos A2/B2/C2 do COVID-19 HG r7, leave_23andme), e
#     (B) catalogo de STRs da coorte (STRling) -> $CATALOG.
#
#   Para cada gene sugestivo, lista os LOCI de STR da coorte cujo gene_name
#   anotado no catalogo coincide com o nome do gene sugestivo. NAO usa
#   sobreposicao de coordenadas — confia na anotacao gene_name do catalogo.
#
#   Gera tres arquivos:
#     1) --out           Dataset unificado: todas as linhas do catalogo + colunas
#                        GWAS (gwas_hit, gwas_p, gwas_phenotypes, gwas_lead_snp).
#     2) --out-overlap   Pares gene x STR (filtrados por gwas_hit=1), para uso
#                        com summarize_subset.py.
#     3) --ids-out       Lista de STRs_ID unicos com hit GWAS (um por linha).
#
# ENTRADAS
#   --catalog        Catalogo de STR da coorte (STRs_analysis_dataset.tsv).
#   --covid-genes    Saida de extract_covid_genes.py filtrada a p<1e-5
#                    (results/covid_genes_suggestive.tsv): chrom, gene_start,
#                    gene_end, gene, best_p, phenotypes, lead_snp.
#   --out            Dataset unificado com colunas GWAS.
#   --out-overlap    Tabela gene x STR para summarize_subset.py.
#   --ids-out        Lista de STRs_ID unicos (um por linha, deduplicado).
#   --debug          Ativa mensagens de debug verbosas em stderr.
# ---------------------------------------------------------------------------
import argparse, csv, sys, time
from collections import Counter


def dbg(msg, debug=False):
    if debug:
        sys.stderr.write(f"[DEBUG] {msg}\n")


def main():
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--covid-genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--out-overlap', required=True)
    ap.add_argument('--ids-out', required=True)
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()
    D = args.debug

    dbg(f"args: catalog={args.catalog} covid-genes={args.covid_genes} "
        f"out={args.out} out-overlap={args.out_overlap} ids-out={args.ids_out}", D)

    # ----------------------------------------------------------------------
    # 1) Carregar genes COVID sugestivos -> dict por nome do gene
    #    genes: {gene_name: [best_p, phenotypes, lead_snp, chrom, gene_start, gene_end]}
    # ----------------------------------------------------------------------
    genes = {}
    with open(args.covid_genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        n_genes_raw = 0
        for row in r:
            n_genes_raw += 1
            gene = row.get('gene', '').strip()
            if not gene:
                continue
            try:
                gs = int(row['gene_start'])
                ge = int(row['gene_end'])
            except (ValueError, KeyError):
                gs, ge = None, None
            chrom = row.get('chrom', '')
            best_p = row.get('best_p', '')
            phenotypes = row.get('phenotypes', '')
            lead_snp = row.get('lead_snp', '')
            genes[gene] = [best_p, phenotypes, lead_snp, chrom, gs, ge]

    sys.stderr.write(f"Genes sugestivos: {len(genes)} unicos de {n_genes_raw} linhas\n")
    if D:
        for g, v in sorted(genes.items()):
            dbg(f"  {g}: p={v[0]}, pheno={v[1]}", D)

    # ----------------------------------------------------------------------
    # 2) Ler catalogo e detectar colunas
    # ----------------------------------------------------------------------
    with open(args.catalog) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    col = {h: i for i, h in enumerate(header)}

    required = ['STRs_ID', 'gene_name']
    missing = [c for c in required if c not in col]
    if missing:
        sys.exit(f"ERRO: colunas nao encontradas no catalogo: {missing}")

    chrom_col = 'chrom' if 'chrom' in col else None
    start_col = 'start' if 'start' in col else None
    end_col = 'end' if 'end' in col else None
    repeat_col = 'repeat_unit' if 'repeat_unit' in col else (
        'repeatunit' if 'repeatunit' in col else None)

    dbg(f"Colunas detectadas: {col}", D)

    # ----------------------------------------------------------------------
    # 3) Primeiro passe: contar carriers por (gene, strs_id)
    #    (cada linha do catalogo = 1 sample carrieira)
    # ----------------------------------------------------------------------
    pair_counts = Counter()
    gene_set = set(genes.keys())

    with open(args.catalog) as fh:
        reader = csv.reader(fh, delimiter='\t')
        next(reader)  # pula header
        for row in reader:
            try:
                gene_name = row[col['gene_name']]
                strs_id = row[col['STRs_ID']]
            except IndexError:
                continue
            if gene_name in gene_set:
                pair_counts[(gene_name, strs_id)] += 1

    sys.stderr.write(f"Pares gene x STR com hit: {len(pair_counts)}\n")

    # ----------------------------------------------------------------------
    # 4) Segundo passe: escrever os 3 arquivos de saida
    # ----------------------------------------------------------------------
    n_total = 0
    n_hit = 0
    idset = set()

    with open(args.catalog) as fh, \
         open(args.out, 'w', newline='') as f_unified, \
         open(args.out_overlap, 'w', newline='') as f_overlap, \
         open(args.ids_out, 'w') as f_ids:

        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)

        # Dataset unificado: header original + 4 colunas GWAS
        w_unified = csv.writer(f_unified, delimiter='\t')
        w_unified.writerow(header + ['gwas_hit', 'gwas_p', 'gwas_phenotypes',
                                      'gwas_lead_snp'])

        # Overlap: gene x STR
        w_overlap = csv.writer(f_overlap, delimiter='\t')
        w_overlap.writerow(['gene', 'chrom', 'gene_start', 'gene_end',
                            'best_p', 'phenotypes', 'lead_snp',
                            'strs_id', 'str_chrom', 'str_start', 'str_end',
                            'repeatunit', 'n_carriers'])

        # Para evitar duplicacao no overlap (1 linha por gene, nao por sample)
        written_pairs = set()

        for row in reader:
            n_total += 1
            strs_id = row[col['STRs_ID']]
            gene_name = row[col['gene_name']]

            if gene_name in gene_set:
                n_hit += 1
                g = genes[gene_name]
                best_p, phenotypes, lead_snp = g[0], g[1], g[2]

                # Linha unificada com colunas GWAS
                w_unified.writerow(row + ['1', best_p, phenotypes, lead_snp])

                # Overlap: 1 linha por (gene, strs_id)
                pair_key = (gene_name, strs_id)
                if pair_key not in written_pairs:
                    written_pairs.add(pair_key)
                    cnt = pair_counts[pair_key]
                    g_chrom = g[3]
                    g_start = g[4] if g[4] is not None else ''
                    g_end = g[5] if g[5] is not None else ''
                    str_chrom = row[col[chrom_col]] if chrom_col else ''
                    str_start = row[col[start_col]] if start_col else ''
                    str_end = row[col[end_col]] if end_col else ''
                    str_repeat = row[col[repeat_col]] if repeat_col else ''

                    w_overlap.writerow([gene_name, g_chrom, g_start, g_end,
                                        best_p, phenotypes, lead_snp,
                                        strs_id, str_chrom, str_start, str_end,
                                        str_repeat, cnt])

                idset.add(strs_id)
                dbg(f"HIT: gene={gene_name} strs_id={strs_id}", D)
            else:
                # Sem hit GWAS — preenche 0 / vazio
                w_unified.writerow(row + ['0', '', '', ''])

    # ----------------------------------------------------------------------
    # 5) Escrever IDs (deduplicado, ja esta no set)
    # ----------------------------------------------------------------------
    with open(args.ids_out, 'w') as fid:
        for sid in sorted(idset):
            fid.write(sid + "\n")

    # ----------------------------------------------------------------------
    # 6) Resumo
    # ----------------------------------------------------------------------
    sys.stderr.write(f"Catalogo: {n_total} linhas totais\n")
    sys.stderr.write(f"STRs em genes sugestivos: {n_hit} linhas ({len(idset)} STRs_ID unicos)\n")
    sys.stderr.write(f"Genes sugestivos com STR na coorte: "
                     f"{len(set(g for g, _ in pair_counts))}/{len(genes)}\n")
    sys.stderr.write(f"Dataset unificado -> {args.out}\n")
    sys.stderr.write(f"Pares gene x STR -> {args.out_overlap}\n")
    sys.stderr.write(f"IDs -> {args.ids_out}\n")
    sys.stderr.write(f"Tempo total: {time.time() - t0:.1f}s\n")


if __name__ == '__main__':
    main()
