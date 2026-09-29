#!/usr/bin/env python3
"""
cross_intervention_STRs.py
--------------------------
Crosses the DEGs of each comparison (TSV file per comparison) with the cohort
STR catalog, including global DBSCAN metrics.

Each GSE subfolder contains one or more DEG TSVs. The comparison name is
extracted from the file name (the DEG(s)_ prefix and suffixes are removed).

Each output row represents ONE STR annotated in ONE DEG gene, for ONE
unique (GSE, comparison) combination. The GSE is kept as a column, so the
same STR can appear N times if its gene is a DEG in N distinct
studies/comparisons.

Significance filters (applied to DEGs):
  --require {fdr,pval,both,either}   what to require (default: fdr)
  --fdr <float>                      FDR threshold         (default: 0.05)
  --pval <float>                     raw p-value threshold (default: 0.05)
  --logfc <float>                    |logFC| threshold     (default: 1.0;
                                     0 disables the filter)

Outputs:
  1) intervention_strs.tsv        - all STRs annotated in DEG genes
                                    (1 row per STR x GSE x comparison)
  2) intervention_outliers.tsv    - only STRs with a global DBSCAN outlier set
  3) intervention_summary.tsv     - summary per (GSE, comparison, gene)
  4) debug_report.txt             - (with --debug) post-execution check

Usage:
  python3 cross_intervention_STRs.py \
      --deg-dir <directory_with_GSE_subfolders> \
      --str-catalog <path/STRs_analysis_dataset.tsv> \
      --out-dir <output_directory> \
      [--require fdr|pval|both|either] \
      [--fdr 0.05] [--pval 0.05] [--logfc 1.0] \
      [--debug]
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# Candidate columns (grouped by type)
# ---------------------------------------------------------------------------
GENE_COL_CANDIDATES  = ['gene_symbol', 'Gene', 'gene', 'gene_name']
FDR_COL_CANDIDATES   = ['FDR', 'adj.P.Val', 'padj', 'p_adj', 'p.adjust']
PVAL_COL_CANDIDATES  = ['P.Value', 'pvalue', 'p_value', 'pval',
                        'PValue', 'p.value']
LOGFC_COL_CANDIDATES = ['logFC', 'logfc', 'log2fc', 'log2FC']
SIG_COL_CANDIDATES   = ['Significant', 'significance']
DIR_COL_CANDIDATES   = ['Direction', 'direction']

SIG_TRUE_VALUES = ('yes', 'true', '1', 'up', 'down',
                   'up_covid', 'down_covid',
                   'up_hfd0', 'down_hfd0',
                   'up_icuvent', 'down_icuvent')


def detect_col(header, candidates):
    low = [h.strip().lower() for h in header]
    for cand in candidates:
        if cand.lower() in low:
            return header[low.index(cand.lower())]
    return None


def extract_intervention(fname):
    """Extract the comparison name from the TSV file name."""
    base = os.path.splitext(fname)[0]
    base = re.sub(r'^GSE\d+_', '', base)
    base = re.sub(r'^DEGs?_', '', base)
    base = re.sub(r'_FDR0\.05_log2FC1$', '', base)
    base = re.sub(r'_ajustado$', '', base)
    base = re.sub(r'_FINAL$', '', base)
    base = re.sub(r'_FDR$', '', base)
    base = re.sub(r'_significativos$', '', base)
    base = re.sub(r'_significativo$', '', base)
    return base


def _to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load_deg_file(path, fdr_thr=0.05, pval_thr=0.05,
                  logfc_thr=1.0, require='fdr', stats=None):
    """
    Read a DEG TSV and return a filtered dict {gene: {...}}.

    require: 'fdr'    -> require FDR < fdr_thr (if FDR is available)
             'pval'   -> require raw p-value < pval_thr
             'both'   -> require both
             'either' -> require at least one of the two
    """
    genes = {}
    if stats is None:
        stats = {}
    for k in ('total_rows', 'sem_gene', 'sig_fail', 'pval_fail',
              'fdr_fail', 'logfc_fail', 'passed'):
        stats.setdefault(k, 0)

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        header = reader.fieldnames
        if not header:
            return genes

        gene_col  = detect_col(header, GENE_COL_CANDIDATES)
        fdr_col   = detect_col(header, FDR_COL_CANDIDATES)
        pval_col  = detect_col(header, PVAL_COL_CANDIDATES)
        logfc_col = detect_col(header, LOGFC_COL_CANDIDATES)
        dir_col   = detect_col(header, DIR_COL_CANDIDATES)
        sig_col   = detect_col(header, SIG_COL_CANDIDATES)

        stats['col_gene']  = gene_col
        stats['col_fdr']   = fdr_col
        stats['col_pval']  = pval_col
        stats['col_logfc'] = logfc_col
        stats['col_sig']   = sig_col

        if gene_col is None:
            sys.stderr.write(
                f"  WARNING: gene column not found in {path}\n")
            return genes

        if require in ('pval', 'both') and pval_col is None:
            sys.stderr.write(
                f"  WARNING: --require {require} but no p-value column in "
                f"{os.path.basename(path)}. p filter ignored.\n")
        if require in ('fdr', 'both') and fdr_col is None:
            sys.stderr.write(
                f"  WARNING: --require {require} but no FDR column in "
                f"{os.path.basename(path)}. FDR filter ignored.\n")

        for row in reader:
            stats['total_rows'] += 1
            gene = row.get(gene_col, '').strip()
            if not gene:
                stats['sem_gene'] += 1
                continue

            # --- filter 1: Significant flag (if present) ---
            if sig_col:
                val = row.get(sig_col, '').strip().lower()
                if val not in SIG_TRUE_VALUES:
                    stats['sig_fail'] += 1
                    continue

            # --- filter 2: raw p-value ---
            pval_ok = True
            if pval_col and require in ('pval', 'both', 'either'):
                p = _to_float(row.get(pval_col, '1'))
                pval_ok = (p is not None) and (p < pval_thr)

            # --- filter 3: FDR ---
            fdr_ok = True
            if fdr_col and require in ('fdr', 'both', 'either'):
                q = _to_float(row.get(fdr_col, '1'))
                fdr_ok = (q is not None) and (q < fdr_thr)

            # --- combination of the two ---
            if require == 'both':
                passed_sig = pval_ok and fdr_ok
            elif require == 'either':
                passed_sig = pval_ok or fdr_ok
            elif require == 'pval':
                passed_sig = pval_ok
            else:  # 'fdr'
                passed_sig = fdr_ok

            if not passed_sig:
                if not pval_ok:
                    stats['pval_fail'] += 1
                if not fdr_ok:
                    stats['fdr_fail'] += 1
                continue

            # --- filter 4: |logFC| ---
            if logfc_col and logfc_thr is not None:
                lfc = _to_float(row.get(logfc_col, '0'))
                if lfc is None or abs(lfc) < logfc_thr:
                    stats['logfc_fail'] += 1
                    continue

            genes[gene] = {
                'logFC':     row.get(logfc_col, '') if logfc_col else '',
                'FDR':       row.get(fdr_col, '')   if fdr_col else '',
                'P.Value':   row.get(pval_col, '')  if pval_col else '',
                'Direction': row.get(dir_col, '')   if dir_col else '',
            }
            stats['passed'] += 1

    return genes


# ---------------------------------------------------------------------------
# DEBUG MODULE
# ---------------------------------------------------------------------------
class DebugReport:
    def __init__(self, enabled):
        self.enabled = enabled
        self.lines = []

    def log(self, msg=''):
        if self.enabled:
            sys.stderr.write(msg + '\n')
        self.lines.append(msg)

    def save(self, path):
        with open(path, 'w') as f:
            f.write('\n'.join(self.lines) + '\n')


def verify_outputs(out_dir, all_matches, outlier_matches, str_catalog,
                   dataset_info, intervention_seen, filter_stats, dbg):
    """Re-read the generated TSVs and check consistency."""
    dbg.log('\n' + '=' * 72)
    dbg.log('POST-EXECUTION CHECK')
    dbg.log('=' * 72)

    # ------------------------------------------------------------------
    # [1] Re-read the generated TSVs
    # ------------------------------------------------------------------
    def read_tsv(path):
        with open(path) as fh:
            return list(csv.DictReader(fh, delimiter='\t'))

    try:
        tsv_strs = read_tsv(os.path.join(out_dir, 'intervention_strs.tsv'))
        tsv_out  = read_tsv(os.path.join(out_dir, 'intervention_outliers.tsv'))
        tsv_sum  = read_tsv(os.path.join(out_dir, 'intervention_summary.tsv'))
    except FileNotFoundError as e:
        dbg.log(f'[ERROR] Output file not found: {e}')
        return

    dbg.log('\n[1] Row counts in the TSVs')
    dbg.log(f'    intervention_strs.tsv      : {len(tsv_strs)}')
    dbg.log(f'    intervention_outliers.tsv  : {len(tsv_out)}')
    dbg.log(f'    intervention_summary.tsv   : {len(tsv_sum)}')
    dbg.log(f'    (in-memory) all_matches    : {len(all_matches)}')
    dbg.log(f'    (in-memory) outlier_matches: {len(outlier_matches)}')

    # ------------------------------------------------------------------
    # [2] In-memory x TSV consistency
    # ------------------------------------------------------------------
    dbg.log('\n[2] In-memory x TSV consistency')
    ok = True
    if len(tsv_strs) != len(all_matches):
        dbg.log(f'    x intervention_strs.tsv has {len(tsv_strs)} rows, '
                f'expected {len(all_matches)}')
        ok = False
    if len(tsv_out) != len(outlier_matches):
        dbg.log(f'    x intervention_outliers.tsv has {len(tsv_out)} rows, '
                f'expected {len(outlier_matches)}')
        ok = False
    if ok:
        dbg.log('    o counts match')

    # ------------------------------------------------------------------
    # [3] GSE-dependence: same STR in same gene across different GSEs
    # ------------------------------------------------------------------
    dbg.log('\n[3] GSE-dependence (STR x gene x GSE)')
    str_gene_gses = defaultdict(set)
    for m in all_matches:
        str_gene_gses[(m['STRs_ID'], m['gene_name'])].add(m['gse'])

    multi_gse = {k: v for k, v in str_gene_gses.items() if len(v) > 1}
    dbg.log(f'    STRs in DEG genes in >=2 GSEs: {len(multi_gse)}')
    if multi_gse:
        dbg.log('    Examples (STR | gene | GSEs):')
        for (sid, gene), gses in sorted(multi_gse.items())[:5]:
            dbg.log(f'      {sid} | {gene} | {sorted(gses)}')
        dbg.log('    Checking that the TSV shows them N times:')
        for (sid, gene), gses in sorted(multi_gse.items())[:3]:
            n_tsv = sum(1 for r in tsv_strs
                        if r['STRs_ID'] == sid and r['gene_name'] == gene)
            status = 'o' if n_tsv == len(gses) else 'x'
            dbg.log(f'      {status} {sid}/{gene}: '
                    f'{n_tsv} rows in the TSV, {len(gses)} GSEs expected')
    else:
        dbg.log('    WARNING: no STR appears in >=2 GSEs - '
                'check whether DEG genes overlap between studies.')

    # ------------------------------------------------------------------
    # [4] Subset: outliers subset of all
    # ------------------------------------------------------------------
    dbg.log('\n[4] outlier_matches subset of all_matches')
    all_keys = Counter((m['gse'], m['intervention'], m['STRs_ID'])
                       for m in all_matches)
    out_keys = Counter((m['gse'], m['intervention'], m['STRs_ID'])
                       for m in outlier_matches)
    inconsistent = 0
    for k, n in out_keys.items():
        if k not in all_keys or n > all_keys[k]:
            inconsistent += 1
    if inconsistent == 0:
        dbg.log(f'    o All {len(outlier_matches)} outliers '
                f'are in all_matches')
    else:
        dbg.log(f'    x {inconsistent} inconsistent outlier keys')

    # ------------------------------------------------------------------
    # [5] Summary uniqueness
    # ------------------------------------------------------------------
    dbg.log('\n[5] (gse, intervention, gene) uniqueness in the summary')
    keys = [(r['gse'], r['intervention'], r['gene']) for r in tsv_sum]
    dups = [k for k, n in Counter(keys).items() if n > 1]
    if not dups:
        dbg.log(f'    o {len(keys)} unique (gse, intervention, gene) pairs')
    else:
        dbg.log(f'    x {len(dups)} duplicated pairs. Examples:')
        for k in dups[:5]:
            dbg.log(f'      {k}')

    # ------------------------------------------------------------------
    # [6] STR catalog coverage
    # ------------------------------------------------------------------
    dbg.log('\n[6] STR catalog coverage')
    catalog_ids  = {r.get('STRs_ID', '').strip() for r in str_catalog}
    annotated_ids = {m['STRs_ID'] for m in all_matches}
    dbg.log(f'    STRs in catalog      : {len(catalog_ids)}')
    dbg.log(f'    STRs annotated (DEGs) : {len(annotated_ids)}')
    dbg.log(f'    STRs without DEG      : {len(catalog_ids - annotated_ids)}')

    # ------------------------------------------------------------------
    # [7] Field integrity
    # ------------------------------------------------------------------
    dbg.log('\n[7] Field integrity in intervention_strs.tsv')
    issues = Counter()
    for r in tsv_strs:
        if not r.get('STRs_ID', '').strip():
            issues['empty STRs_ID'] += 1
        if not r.get('gene_name', '').strip():
            issues['empty gene_name'] += 1
        if not r.get('gse', '').strip():
            issues['empty gse'] += 1
        if not r.get('intervention', '').strip():
            issues['empty intervention'] += 1
        lfc = r.get('logFC', '').strip()
        if lfc and _to_float(lfc) is None:
            issues['non-numeric logFC'] += 1
        fdr = r.get('FDR', '').strip()
        if fdr and _to_float(fdr) is None:
            issues['non-numeric FDR'] += 1
    if not issues:
        dbg.log('    o No integrity problems')
    else:
        for k, v in issues.items():
            dbg.log(f'    x {k}: {v} rows')

    # ------------------------------------------------------------------
    # [8] intervention collisions between GSEs
    # ------------------------------------------------------------------
    dbg.log('\n[8] intervention collisions between GSEs')
    collisions = {k: v for k, v in intervention_seen.items() if len(v) > 1}
    if not collisions:
        dbg.log('    o No intervention shared between GSEs')
    else:
        dbg.log(f'    WARNING: {len(collisions)} interventions in multiple GSEs:')
        for interv, gses in sorted(collisions.items()):
            dbg.log(f'      {interv}: {sorted(gses)}')

    # ------------------------------------------------------------------
    # [9] Distribution by (gse, intervention)
    # ------------------------------------------------------------------
    dbg.log('\n[9] Row distribution by (gse, intervention)')
    per_ds = Counter((m['gse'], m['intervention']) for m in all_matches)
    for (gse, interv), n in sorted(per_ds.items()):
        dbg.log(f'    {gse} | {interv}: {n} STR rows')

    # ------------------------------------------------------------------
    # [10] Genes per (gse, intervention) - vs dataset_info
    # ------------------------------------------------------------------
    dbg.log('\n[10] DEG genes per (gse, intervention) - vs dataset_info')
    per_ds_genes = defaultdict(set)
    for m in all_matches:
        per_ds_genes[(m['gse'], m['intervention'])].add(m['gene_name'])
    for key, genes in sorted(per_ds_genes.items()):
        info = dataset_info.get(key, {})
        n_expected = info.get('n_degs', '?')
        dbg.log(f'    {key[0]} | {key[1]}: {len(genes)} genes with STR '
                f'(filtered DEGs in file: {n_expected})')

    # ------------------------------------------------------------------
    # [11] Spot-check: does logFC vary by GSE?
    # ------------------------------------------------------------------
    dbg.log('\n[11] Spot-check: does logFC differ between GSEs for the same gene?')
    gene_vals = defaultdict(dict)
    for m in all_matches:
        gene_vals[m['gene_name']][m['gse']] = m['logFC']
    differing = 0
    same = 0
    for gene, vals in gene_vals.items():
        if len(vals) > 1:
            unique_vals = {v for v in vals.values() if v}
            if len(unique_vals) > 1:
                differing += 1
            else:
                same += 1
    dbg.log(f'    Genes with different logFC between GSEs: {differing}')
    dbg.log(f'    Genes with identical logFC between GSEs : {same}')
    if differing == 0 and same > 0:
        dbg.log('    WARNING: no gene had a distinct logFC between GSEs - '
                'check whether values are being overwritten.')
    elif differing > 0:
        dbg.log('    o logFC varies by GSE (GSE-dependence confirmed)')

    # ------------------------------------------------------------------
    # [12] DEG filtering statistics
    # ------------------------------------------------------------------
    dbg.log('\n[12] DEG filtering by (gse, intervention)')
    if not filter_stats:
        dbg.log('    (no filtering data)')
    else:
        dbg.log('    dataset | rows | sig_fail | pval_fail | fdr_fail | '
                'logfc_fail | passed | col_pval | col_fdr')
        for key in sorted(filter_stats):
            s = filter_stats[key]
            dbg.log(
                f"    {key[0]} | {key[1]} | "
                f"rows={s.get('total_rows', 0)} | "
                f"sig_fail={s.get('sig_fail', 0)} | "
                f"pval_fail={s.get('pval_fail', 0)} | "
                f"fdr_fail={s.get('fdr_fail', 0)} | "
                f"logfc_fail={s.get('logfc_fail', 0)} | "
                f"passed={s.get('passed', 0)} | "
                f"col_pval={s.get('col_pval')} | "
                f"col_fdr={s.get('col_fdr')}"
            )

    dbg.log('\n' + '=' * 72)
    dbg.log('END OF CHECK')
    dbg.log('=' * 72)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description='Crosses DEGs per comparison with cohort STRs')
    ap.add_argument('--deg-dir', required=True,
                    help='Root directory with GSE subfolders')
    ap.add_argument('--str-catalog', required=True,
                    help='Path to STRs_analysis_dataset.tsv')
    ap.add_argument('--out-dir', default='.',
                    help='Output directory')
    ap.add_argument('--fdr', type=float, default=0.05,
                    help='FDR threshold (default: 0.05)')
    ap.add_argument('--pval', type=float, default=0.05,
                    help='Raw p-value threshold (default: 0.05)')
    ap.add_argument('--logfc', type=float, default=1.0,
                    help='|logFC| threshold (default: 1.0). '
                         'Use 0 to disable.')
    ap.add_argument('--require', choices=['fdr', 'pval', 'both', 'either'],
                    default='fdr',
                    help='Which significance tests to require (default: fdr)')
    ap.add_argument('--debug', action='store_true',
                    help='Run post-execution check and generate debug_report.txt')
    args = ap.parse_args()

    dbg = DebugReport(enabled=args.debug)
    os.makedirs(args.out_dir, exist_ok=True)

    logfc_thr = args.logfc if args.logfc > 0 else None
    sys.stderr.write(
        f"DEG filters: require={args.require}, "
        f"pval<{args.pval}, fdr<{args.fdr}, "
        f"|logFC|>{logfc_thr if logfc_thr is not None else 'off'}\n")

    gse_dirs = sorted([d for d in glob.glob(os.path.join(args.deg_dir, 'GSE*'))
                        if os.path.isdir(d)])
    if not gse_dirs:
        sys.stderr.write(f"ERROR: no GSE folder in {args.deg_dir}\n")
        sys.exit(1)
    sys.stderr.write(f"GSE folders found: {len(gse_dirs)}\n")

    sys.stderr.write(f"Loading catalog: {args.str_catalog}\n")
    str_catalog = []
    with open(args.str_catalog) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            str_catalog.append(row)
    sys.stderr.write(f"  {len(str_catalog)} STRs loaded\n")

    all_matches = []
    outlier_matches = []
    dataset_info = {}
    intervention_seen = {}
    filter_stats = {}

    for gse_dir in gse_dirs:
        gse_name = os.path.basename(gse_dir)
        deg_files = sorted(glob.glob(os.path.join(gse_dir, '*.tsv')))

        for deg_file in deg_files:
            fname = os.path.basename(deg_file)
            intervention = extract_intervention(fname)
            dataset_label = f"{gse_name}/{fname}"
            key = (gse_name, intervention)

            sys.stderr.write(f"\n=== {dataset_label} "
                             f"(comparison: {intervention}) ===\n")

            degs_stats = {}
            degs = load_deg_file(
                deg_file,
                fdr_thr=args.fdr,
                pval_thr=args.pval,
                logfc_thr=logfc_thr,
                require=args.require,
                stats=degs_stats,
            )
            filter_stats[key] = degs_stats

            sys.stderr.write(
                f"  {len(degs)} DEG genes after filter "
                f"(rows={degs_stats.get('total_rows', 0)}, "
                f"sig_fail={degs_stats.get('sig_fail', 0)}, "
                f"pval_fail={degs_stats.get('pval_fail', 0)}, "
                f"fdr_fail={degs_stats.get('fdr_fail', 0)}, "
                f"logfc_fail={degs_stats.get('logfc_fail', 0)})\n")

            if not degs:
                continue

            intervention_seen.setdefault(intervention, set()).add(gse_name)
            if len(intervention_seen[intervention]) > 1:
                sys.stderr.write(
                    f"  WARNING: comparison '{intervention}' in multiple GSEs: "
                    f"{sorted(intervention_seen[intervention])}\n")

            dataset_info[key] = {
                'gse': gse_name,
                'file': fname,
                'intervention': intervention,
                'n_degs': len(degs),
            }

            n_all = 0
            n_outlier = 0
            for str_row in str_catalog:
                gene_name = str_row.get('gene_name', '').strip()
                if gene_name not in degs:
                    continue
                sid = str_row.get('STRs_ID', '').strip()
                match = {
                    'gse': gse_name,
                    'intervention': intervention,
                    'dataset': dataset_label,
                    'gene_name': gene_name,
                    'STRs_ID': sid,
                    'chrom': str_row.get('chrom', ''),
                    'start': str_row.get('start', ''),
                    'end': str_row.get('end', ''),
                    'repeat_unit': str_row.get('repeat_unit', ''),
                    'allele1_est': str_row.get('allele1_est', ''),
                    'allele2_est': str_row.get('allele2_est', ''),
                    'depth': str_row.get('depth', ''),
                    'region': str_row.get('region', ''),
                    'group': str_row.get('group', ''),
                    'logFC': degs[gene_name]['logFC'],
                    'FDR': degs[gene_name]['FDR'],
                    'P.Value': degs[gene_name]['P.Value'],
                    'Direction': degs[gene_name]['Direction'],
                    'n_outliers_dbscan_global':
                        str_row.get('n_outliers_dbscan_global', ''),
                    'outlier_samples_dbscan_global':
                        str_row.get('outlier_samples_dbscan_global', ''),
                    'n_clusters_dbscan_global':
                        str_row.get('n_clusters_dbscan_global', ''),
                    'noise_ratio_dbscan_global':
                        str_row.get('noise_ratio_dbscan_global', ''),
                }
                all_matches.append(match)
                n_all += 1

                try:
                    n_clusters = int(match['n_clusters_dbscan_global'])
                    noise_ratio = float(match['noise_ratio_dbscan_global'])
                    n_outliers_val = int(match['n_outliers_dbscan_global'])
                except (ValueError, KeyError):
                    continue
                if n_clusters > 0 and noise_ratio <= 0.10 and n_outliers_val >= 1:
                    outlier_matches.append(match)
                    n_outlier += 1

            sys.stderr.write(f"  STRs: {n_all} total, {n_outlier} outliers\n")

    # ----------------------------------------------------------------------
    # OUTPUTS
    # ----------------------------------------------------------------------
    fields = ['gse', 'intervention', 'dataset', 'gene_name', 'STRs_ID',
              'chrom', 'start', 'end', 'repeat_unit',
              'allele1_est', 'allele2_est', 'depth', 'region', 'group',
              'logFC', 'FDR', 'P.Value', 'Direction',
              'n_outliers_dbscan_global', 'outlier_samples_dbscan_global',
              'n_clusters_dbscan_global', 'noise_ratio_dbscan_global']

    out_strs = os.path.join(args.out_dir, 'intervention_strs.tsv')
    with open(out_strs, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(all_matches)
    sys.stderr.write(f"\nWritten: {out_strs} ({len(all_matches)} rows)\n")

    out_out = os.path.join(args.out_dir, 'intervention_outliers.tsv')
    with open(out_out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(outlier_matches)
    sys.stderr.write(f"Written: {out_out} ({len(outlier_matches)} rows)\n")

    def largest_allele(m):
        a1 = _to_float(m.get('allele1_est'))
        a2 = _to_float(m.get('allele2_est'))
        if a1 is None and a2 is None:
            return None
        if a1 is None:
            return a2
        if a2 is None:
            return a1
        return max(a1, a2)

    locus_group_vals = {}
    for m in all_matches:
        sid = m['STRs_ID']
        g = (m.get('group') or '').strip()
        big = largest_allele(m)
        if not g or big is None:
            continue
        locus_group_vals.setdefault(sid, {}).setdefault(g, []).append(big)

    def locus_overlap(sid):
        groups = locus_group_vals.get(sid, {})
        labels = [g for g, vals in groups.items() if vals]
        if len(labels) < 2:
            return None
        ranges = {g: (min(groups[g]), max(groups[g])) for g in labels}
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a = ranges[labels[i]]
                b = ranges[labels[j]]
                if not (a[0] <= b[1] and b[0] <= a[1]):
                    return False
        return True

    gene_str_loci = {}
    for m in all_matches:
        gene_str_loci.setdefault(m['gene_name'], set()).add(m['STRs_ID'])

    def gene_overlap(gene):
        flags = [locus_overlap(sid) for sid in gene_str_loci.get(gene, ())]
        flags = [f for f in flags if f is not None]
        if not flags:
            return 'sem_dados'
        return 'nao' if any(f is False for f in flags) else 'sim'

    intervention_keys = sorted(set(
        (m['gse'], m['intervention']) for m in all_matches))

    outlier_gene_pairs = {}
    for m in outlier_matches:
        k = (m['gse'], m['intervention'])
        outlier_gene_pairs.setdefault(k, set()).add(
            (m['gene_name'], m['STRs_ID']))

    summary_rows = []
    for (gse_name, interv) in intervention_keys:
        interv_matches = [m for m in all_matches
                          if m['gse'] == gse_name
                          and m['intervention'] == interv]
        gene_loci = {}
        for m in interv_matches:
            gene_loci.setdefault(m['gene_name'], set()).add(m['STRs_ID'])
        for gene in sorted(gene_loci):
            loci = gene_loci[gene]
            outlier_loci = {sid for (gn, sid)
                            in outlier_gene_pairs.get((gse_name, interv), ())
                            if gn == gene}
            summary_rows.append({
                'gse': gse_name,
                'intervention': interv,
                'gene': gene,
                'n_strs_identified': len(loci),
                'n_strs_identified_outliers': len(outlier_loci),
                'overlap_maior_alealo_grupos': gene_overlap(gene),
            })

    out_sum = os.path.join(args.out_dir, 'intervention_summary.tsv')
    with open(out_sum, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gse', 'intervention', 'gene',
                                          'n_strs_identified',
                                          'n_strs_identified_outliers',
                                          'overlap_maior_alealo_grupos'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)
    sys.stderr.write(f"Written: {out_sum} ({len(summary_rows)} rows)\n")

    sys.stderr.write("\n=== Summary per (GSE, comparison) ===\n")
    for (gse_name, interv) in intervention_keys:
        rows = [r for r in summary_rows
                if r['gse'] == gse_name and r['intervention'] == interv]
        n_strs = sum(r['n_strs_identified'] for r in rows)
        n_out = sum(r['n_strs_identified_outliers'] for r in rows)
        n_nao = sum(1 for r in rows
                    if r['overlap_maior_alealo_grupos'] == 'nao')
        sys.stderr.write(
            f"  {gse_name} | {interv}: {n_strs} STRs ({n_out} outliers), "
            f"{len(rows)} genes; {n_nao} WITHOUT overlap\n")

    # ----------------------------------------------------------------------
    # DEBUG / CHECK
    # ----------------------------------------------------------------------
    if args.debug:
        verify_outputs(
            out_dir=args.out_dir,
            all_matches=all_matches,
            outlier_matches=outlier_matches,
            str_catalog=str_catalog,
            dataset_info=dataset_info,
            intervention_seen=intervention_seen,
            filter_stats=filter_stats,
            dbg=dbg,
        )
        report_path = os.path.join(args.out_dir, 'debug_report.txt')
        dbg.save(report_path)
        sys.stderr.write(f'\nDebug report: {report_path}\n')

    sys.stderr.write("\nDone.\n")


if __name__ == '__main__':
    main()