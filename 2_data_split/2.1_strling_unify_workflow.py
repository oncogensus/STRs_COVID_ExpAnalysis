#!/usr/bin/env python
import argparse
import sys
import os
import pandas as pd
import datetime
import time

def log_debug(message, start_time=None):
    """Prints debug logs with timestamp and optional elapsed time."""
    now = datetime.datetime.now().strftime("%H:%M:%S")
    elapsed = f" [Elapsed: {datetime.timedelta(seconds=int(time.time() - start_time))}]" if start_time else ""
    print(f"[{now}]{elapsed} DEBUG: {message}", file=sys.stderr)

def parse_args():
    parser = argparse.ArgumentParser(description='STRling Workflow: Quality Filtering and Summary Generation')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory with STRling -genotype.txt files')
    parser.add_argument('--groups', type=str, required=True, help='CSV with sample/group info')
    parser.add_argument('--out_dir', type=str, required=True, help='Output directory')
    # Quality parameters
    parser.add_argument('--min_depth', type=float, default=15.0, help='Filter: minimum depth (default: 15)')
    parser.add_argument('--min_clips_sum', type=int, default=1, help='Filter: minimum sum of L+R clips (default > 0)')
    parser.add_argument('--exclude_homopolymers', action='store_true', default=True, help='Exclude loci with repeat unit length of 1')
    return parser.parse_args()

def main():
    start_time = time.time()
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    log_debug(f"Starting workflow: Filtering for Depth >= {args.min_depth}, Clips Sum >= {args.min_clips_sum}, Exclude Homopolymers: {args.exclude_homopolymers}", start_time)

    # 1. Load Group Metadata
    try:
        groups_df = pd.read_csv(args.groups)
        groups_df['sample'] = groups_df['sample'].astype(str).str.strip()
        sample_list = groups_df['sample'].tolist()
    except Exception as e:
        sys.exit(f"ERROR: Could not read groups file: {e}")

    genotype_dfs = []
    
    # 2. Load Genotype Files
    log_debug("Loading genotype files...")
    for sample in sample_list:
        gt_path = os.path.join(args.input_dir, f"{sample}-genotype.txt")
        if os.path.exists(gt_path):
            temp_df = pd.read_csv(gt_path, sep='\t', na_values="na", low_memory=False)
            if "#chrom" in temp_df.columns:
                temp_df.rename(columns={"#chrom": "chrom"}, inplace=True)
            temp_df['sample'] = sample
            genotype_dfs.append(temp_df)
        else:
            log_debug(f"WARNING: File missing for sample {sample}")

    if not genotype_dfs:
        sys.exit("ERROR: No valid genotype files found.")

    # Initial Merge
    full_data = pd.concat(genotype_dfs, axis=0, ignore_index=True, sort=False)
    full_data = full_data.merge(groups_df, on="sample", how="left")

    # 3. Preparation and Numeric Conversion
    numeric_cols = ['left_clips', 'right_clips', 'depth', 'allele1_est', 'allele2_est', 'spanning_reads']
    for col in numeric_cols:
        if col in full_data.columns:
            full_data[col] = pd.to_numeric(full_data[col], errors='coerce').fillna(0)

    # Create Locus ID: chrom-left-right-repeatunit
    full_data['locus'] = (
        full_data['chrom'].astype(str) + "-" + 
        full_data['left'].astype(str).str.replace(r'\.0$', '', regex=True) + "-" + 
        full_data['right'].astype(str).str.replace(r'\.0$', '', regex=True) + "-" + 
        full_data['repeatunit'].astype(str)
    )

    # 4. Quality Filtering Application
    initial_count = len(full_data)
    
    # Filter 4a: Homopolymers
    if args.exclude_homopolymers:
        # Repeats where the unit length is 1 (A, C, T, G)
        full_data = full_data[full_data['repeatunit'].astype(str).str.len() > 1].copy()
        after_hp_count = len(full_data)
        log_debug(f"Homopolymers removed: {initial_count - after_hp_count} rows.")
    else:
        after_hp_count = initial_count

    # Filter 4b: Depth and Clips
    full_data['total_clips'] = full_data['left_clips'] + full_data['right_clips']
    full_data = full_data[
        (full_data['depth'] >= args.min_depth) & 
        (full_data['total_clips'] >= args.min_clips_sum)
    ].copy()
    
    log_debug(f"Quality filters applied: {after_hp_count} -> {len(full_data)} rows retained.")
    log_debug(f"Total rows filtered out: {initial_count - len(full_data)}")

    # 5. Zygosity Logic
    full_data['allele1_rounded'] = full_data['allele1_est'].round()
    full_data['allele2_rounded'] = full_data['allele2_est'].round()
    full_data['is_homozygous'] = (full_data['allele1_rounded'] == full_data['allele2_rounded']).astype(int)
    full_data['is_heterozygous'] = (full_data['allele1_rounded'] != full_data['allele2_rounded']).astype(int)

    # 6. SUMMARY GENERATION
    log_debug("Generating summary statistics...")
    patient_summary = full_data.groupby(['sample', 'group']).agg(
        variant_count=('locus', 'count'),
        avg_depth=('depth', 'mean'),
        homozygous_count=('is_homozygous', 'sum'),
        heterozygous_count=('is_heterozygous', 'sum')
    ).reset_index()

    group_summary = patient_summary.groupby('group').agg(
        total_variants=('variant_count', 'sum'),
        avg_variants_per_patient=('variant_count', 'mean'),
        avg_depth=('avg_depth', 'mean'),
        total_homozygous=('homozygous_count', 'sum'),
        total_heterozygous=('heterozygous_count', 'sum')
    ).reset_index()

    global_metrics = pd.DataFrame([{
        'group': 'GLOBAL_TOTAL',
        'total_variants': patient_summary['variant_count'].sum(),
        'avg_variants_per_patient': patient_summary['variant_count'].mean(),
        'avg_depth': patient_summary['avg_depth'].mean(),
        'total_homozygous': patient_summary['homozygous_count'].sum(),
        'total_heterozygous': patient_summary['heterozygous_count'].sum()
    }])

    # 7. EXPORTING FILES
    log_debug("Exporting result files...")
    full_data.to_csv(os.path.join(args.out_dir, "global_STRs_filtered.tsv"), sep='\t', index=False)
    patient_summary.to_csv(os.path.join(args.out_dir, "summary_by_patient.tsv"), sep='\t', index=False)
    final_report = pd.concat([group_summary, global_metrics], ignore_index=True)
    final_report.to_csv(os.path.join(args.out_dir, "summary_report_final.tsv"), sep='\t', index=False)

    log_debug(f"Workflow finished successfully. Results in: {args.out_dir}", start_time)

if __name__ == '__main__':
    main()