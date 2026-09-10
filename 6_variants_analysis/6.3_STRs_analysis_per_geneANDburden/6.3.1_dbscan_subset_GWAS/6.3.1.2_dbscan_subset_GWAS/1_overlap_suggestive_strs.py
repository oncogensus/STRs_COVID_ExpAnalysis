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
#   Suporta multiplos hits por gene (multiplas linhas no covid_genes_suggestive.tsv).
#   Para cada gene, agrega:
#     - best_p: menor p-value entre todos os SNPs do gene
#     - all_snps_p: todos os p-values separados por ;
#     - all_snps_ids: todos os rsids separados por ;
#     - phenotypes: uniao de todos os fenotipos
#     - lead_snp: rsid com melhor p-value
#
#   Gera tres arquivos:
#     1) --out           Dataset unificado: todas as linhas do catalogo + colunas
#                        GWAS (gwas_hit, gwas_p, gwas_all_p, gwas_phenotypes,
#                        gwas_lead_snp).
#     2) --out-overlap   Pares gene x STR (filtrados por gwas_hit=1), com
#                        best_p + all_snps_p + all_snps_ids.
#     3) --ids-out       Lista de STRs_ID unicos com hit GWAS (um por linha).
#
# ENTRADAS
#   --catalog        Catalogo de STR da coorte (STRs_analysis_dataset.tsv).
#   --covid-genes    Saida de extract_covid_genes.py filtrada a p<1e-5
#                    (results/covid_genes_suggestive.tsv): gene, chrom, gene_start,
#                    gene_end, snp_p, snp_phenotype, snp_rsid.
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
    # 1) Carregar genes COVID sugestivos
#    genes: {gene_name: [(snp_p, snp_phenotype, snp_rsid, chrom, gene_start, gene_end), ...]}
#    Cada linha do arquivo = 1 SNP hit. Multiplas linhas por gene = multiplos hits.
    # ----------------------------------------------------------------------
    genes = {}  # gene -> lista de (p, pheno, rsid, chrom, gs, ge)
    with open(args.covid_genes) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        n_rows_raw = 0
        for row in r:
            n_rows_raw += 1
            gene = row.get('gene', '').strip()
            if not gene:
                continue
            try:
                p = float(row.get('snp_p', ''))
            except (ValueError, TypeError):
                continue
            pheno = row.get('snp_phenotype', '')
            rsid = row.get('snp_rsid', '')
            chrom = row.get('chrom', '')
            try:
                gs = int(row['gene_start'])
                ge = int(row['gene_end'])
            except (ValueError, KeyError):
                gs, ge = None, None
            genes.setdefault(gene, []).append((p, pheno, rsid, chrom, gs, ge))

    # Agregar informacao por gene
    # gene_info: {gene: [best_p, all_p_str, phenotypes, lead_snp, chrom, gs, ge]}
    gene_info = {}
    for gene, hits in genes.items():
        hits_sorted = sorted(hits, key=lambda x: x[0])
        best_p = hits_sorted[0][0]
        all_p = ";".join(f"{h[0]:.3e}" for h in hits_sorted)
        all_rsids = ";".join(h[2] for h in hits_sorted if h[2])
        phenotypes = ",".join(sorted(set(h[1] for h in hits if h[1])))
        lead_snp = hits_sorted[0][2]  # rsid com melhor p
        chrom = hits_sorted[0][3]
        gs = hits_sorted[0][4]
        ge = hits_sorted[0][5]
        gene_info[gene] = [best_p, all_p, phenotypes, lead_snp, all_rsids, chrom, gs, ge]

    n_genes = len(gene_info)
    n_hits = sum(len(v) for v in genes.values())
    sys.stderr.write(f"Genes sugestivos: {n_genes} unicos de {n_rows_raw} linhas "
                     f"({n_hits} SNP hits total)\n")
    if D:
        for g, v in sorted(gene_info.items()):
            dbg(f"  {g}: best_p={v[0]:.3e}, hits={len(genes[g])}, pheno={v[2]}", D)

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
    # ----------------------------------------------------------------------
    pair_counts = Counter()
    gene_set = set(gene_info.keys())

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

        # Dataset unificado: header original + 5 colunas GWAS
        w_unified = csv.writer(f_unified, delimiter='\t')
        w_unified.writerow(header + ['gwas_hit', 'gwas_p', 'gwas_all_p',
                                      'gwas_phenotypes', 'gwas_lead_snp'])

        # Overlap: gene x STR
        w_overlap = csv.writer(f_overlap, delimiter='\t')
        w_overlap.writerow(['gene', 'chrom', 'gene_start', 'gene_end',
                            'best_p', 'all_snps_p', 'all_snps_ids',
                            'phenotypes', 'lead_snp',
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
                gi = gene_info[gene_name]
                best_p, all_p, phenotypes, lead_snp = gi[0], gi[1], gi[2], gi[3]

                # Linha unificada com colunas GWAS
                w_unified.writerow(row + ['1', f"{best_p:.3e}", all_p,
                                          phenotypes, lead_snp])

                # Overlap: 1 linha por (gene, strs_id)
                pair_key = (gene_name, strs_id)
                if pair_key not in written_pairs:
                    written_pairs.add(pair_key)
                    cnt = pair_counts[pair_key]
                    g_chrom = gi[5]
                    g_start = gi[6] if gi[6] is not None else ''
                    g_end = gi[7] if gi[7] is not None else ''
                    str_chrom = row[col[chrom_col]] if chrom_col else ''
                    str_start = row[col[start_col]] if start_col else ''
                    str_end = row[col[end_col]] if end_col else ''
                    str_repeat = row[col[repeat_col]] if repeat_col else ''

                    w_overlap.writerow([gene_name, g_chrom, g_start, g_end,
                                        f"{best_p:.3e}", all_p, gi[4],
                                        phenotypes, lead_snp,
                                        strs_id, str_chrom, str_start, str_end,
                                        str_repeat, cnt])

                idset.add(strs_id)
                dbg(f"HIT: gene={gene_name} strs_id={strs_id}", D)
            else:
                # Sem hit GWAS — preenche 0 / vazio
                w_unified.writerow(row + ['0', '', '', '', ''])

    # ----------------------------------------------------------------------
    # 5) Escrever IDs (deduplicado)
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
                     f"{len(set(g for g, _ in pair_counts))}/{len(gene_info)}\n")
    sys.stderr.write(f"Dataset unificado -> {args.out}\n")
    sys.stderr.write(f"Pares gene x STR -> {args.out_overlap}\n")
    sys.stderr.write(f"IDs -> {args.ids_out}\n")
    sys.stderr.write(f"Tempo total: {time.time() - t0:.1f}s\n")


if __name__ == '__main__':
    main()
