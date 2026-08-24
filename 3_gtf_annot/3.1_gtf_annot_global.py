"""
STR Annotation Pipeline - Python Version
Features:
1. Detailed annotation for protein-coding genes (CDS, UTR, Exons, Introns, Promoters).
2. Non-coding genes classified as 'others'.
3. Empty genomic space classified as 'intergenic'.
"""

import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict

import polars as pl
from pybedtools import BedTool

# ==============================
# CONFIGURATION
# ==============================

class Config:
    """Pipeline settings"""
    TSV_PATH = "../samples/global_STRs_filtered.tsv"
    GTF_PATH = "Homo_sapiens.GRCh38.98.gtf"
    GENOME_FILE = "genome.txt"
    OUTDIR = "../samples/"
    PROMOTER_DIST = 3000
    
    # Priority hierarchy - UPDATED to include 'exon'
    PRIORITY_MAP = {
        "CDS": 1,
        "five_prime_utr": 2,
        "three_prime_utr": 3,
        "non_coding_exons": 4,           # Non-translated exons (or specific isoforms)
        "promoter": 5,
        "intron": 6,
        "others": 7,         # Non-coding genes
        "intergenic": 8      # Genomic regions between genes
    }

    GENE_ID_COLUMN = "gene_id"
    KEEP_TEMP_FILES = False

# Create output directory
Path(Config.OUTDIR).mkdir(parents=True, exist_ok=True)

# ==============================
# HELPER FUNCTIONS
# ==============================

def add_chr(chrom: str) -> str:
    """Add 'chr' to chromosome if not present"""
    chrom_str = str(chrom)
    if chrom_str.startswith("chr"):
        return chrom_str
    return f"chr{chrom_str}"

def validate_input_files():
    """Validate existence of input files"""
    for path, desc in [
        (Config.TSV_PATH, "TSV file"),
        (Config.GTF_PATH, "GTF file"),
        (Config.GENOME_FILE, "Genome file")
    ]:
        if not Path(path).exists():
            raise FileNotFoundError(f"{desc} not found: {path}")

def clean_temp_files():
    """Remove temporary files"""
    if not Config.KEEP_TEMP_FILES:
        temp_files = [
            "coding_genes.bed", "non_coding_genes.bed", "all_genes.bed",
            "exons.bed", "CDS.bed", "five_prime_utr.bed", "three_prime_utr.bed",
            "promoters.bed", "introns.bed", "intergenic.bed", "str_regions_temp.bed"
        ]
        for file in temp_files:
            path = Path(Config.OUTDIR) / file
            if path.exists():
                try: path.unlink()
                except: pass

def read_gtf_features(gtf_path: str) -> pl.DataFrame:
    """Reads GTF and extracts biotype to distinguish coding from non-coding."""
    print(f"Reading GTF: {gtf_path}")
    schema = {
        "chrom": pl.Utf8, "source": pl.Utf8, "feature": pl.Utf8,
        "start": pl.Int64, "end": pl.Int64, "score": pl.Utf8,
        "strand": pl.Utf8, "frame": pl.Utf8, "attribute": pl.Utf8
    }
    df = pl.read_csv(
        gtf_path, separator="\t", has_header=False,
        comment_prefix="#", null_values=".",
        new_columns=list(schema.keys()),
        schema_overrides=schema,
        infer_schema_length=0
    )
    df = df.with_columns([
        pl.col("attribute").str.extract(r'gene_id\s+"([^"]+)"', 1).alias("gene_id"),
        pl.col("attribute").str.extract(r'gene_name\s+"([^"]+)"', 1).alias("gene_name"),
        pl.col("attribute").str.extract(r'gene_biotype\s+"([^"]+)"', 1).alias("gene_biotype")
    ])
    df = df.with_columns(
        pl.when(pl.col("gene_name").is_null())
        .then(pl.col("gene_id"))
        .otherwise(pl.col("gene_name"))
        .alias("gene_name")
    ).with_columns(
        pl.col("chrom").map_elements(add_chr, return_dtype=pl.Utf8).alias("chrom")
    )
    return df

def get_valid_chromosomes(genome_file: str) -> Set[str]:
    """Read genome file and return valid chromosomes"""
    valid_chroms = set()
    with open(genome_file, 'r') as f:
        for line in f:
            if line.strip():
                valid_chroms.add(line.split('\t')[0].strip())
    return valid_chroms

def filter_valid_chromosomes(df: pl.DataFrame, valid_chroms: Set[str]) -> pl.DataFrame:
    """Filter DataFrame to keep only valid chromosomes"""
    return df.filter(pl.col("chrom").is_in(list(valid_chroms)))

def save_bedtool(df: pl.DataFrame, filename: str, name_col: str = None) -> BedTool:
    """Save DataFrame as BED file and return BedTool"""
    path = Path(Config.OUTDIR) / filename
    if name_col is None: name_col = Config.GENE_ID_COLUMN
    bed_df = df.with_columns([
        (pl.col("start") - 1).alias("start"),
        pl.col("end").alias("end")
    ]).select(["chrom", "start", "end", pl.col(name_col).alias("name")])
    bed_df.write_csv(path, separator="\t", include_header=False)
    return BedTool(str(path))

def create_promoters(genes_bed: BedTool) -> BedTool:
    """Create promoter regions upstream of TSS"""
    if genes_bed.count() == 0: return BedTool("", from_string=True)
    tss_records = [[g.chrom, g.start, g.start + 1, g.name] for g in genes_bed]
    tss_bed = BedTool("\n".join("\t".join(map(str, r)) for r in tss_records), from_string=True)
    promoters = tss_bed.slop(g=Config.GENOME_FILE, l=Config.PROMOTER_DIST, r=0, s=False).sort()
    promoters.saveas(str(Path(Config.OUTDIR) / "promoters.bed"))
    return promoters

def create_introns(genes_bed: BedTool, exons_bed: BedTool) -> BedTool:
    """Create introns by subtracting exons from genes"""
    if genes_bed.count() == 0: return BedTool("", from_string=True)
    introns = genes_bed.subtract(exons_bed.sort()).sort()
    introns.saveas(str(Path(Config.OUTDIR) / "introns.bed"))
    return introns

def process_tsv_to_bed(tsv_path: str) -> Tuple[BedTool, pl.DataFrame]:
    """Process consolidated TSV file into BED format"""
    df_tsv = pl.read_csv(tsv_path, separator="\t", has_header=True, infer_schema_length=10000)
    
    # Standardize column names
    rename_map = {"left": "start", "right": "end", "repeatunit": "repeat_unit"}
    df_tsv = df_tsv.rename({k: v for k, v in rename_map.items() if k in df_tsv.columns})

    # Prepare BED-compatible columns
    df = df_tsv.with_row_index("row_id").with_columns([
        pl.col("chrom").map_elements(add_chr, return_dtype=pl.Utf8).alias("var_chrom"),
        (pl.col("start") - 1).alias("var_start"),
        pl.col("end").alias("var_end"),
        pl.col("row_id").cast(pl.Utf8).alias("key"),
        (pl.col("end") - pl.col("start") + 1).alias("region_length"),
        pl.col("sample").alias("sample_id") if "sample" in df_tsv.columns else pl.lit("unknown").alias("sample_id")
    ])
    
    df = df.with_columns([
        (pl.col("region_length") / pl.col("repeat_unit").str.len_bytes()).alias("repeat_count")
    ]).with_columns([
        (pl.col("var_chrom") + ":" + pl.col("start").cast(pl.Utf8) + ":" + 
         pl.col("repeat_unit") + ":" + pl.col("repeat_count").cast(pl.Int64).cast(pl.Utf8)).alias("STRs_ID")
    ])

    bed_path = Path(Config.OUTDIR) / "str_regions_temp.bed"
    df.select(["var_chrom", "var_start", "var_end", "key"]).write_csv(bed_path, separator="\t", include_header=False)
    return BedTool(str(bed_path)).sort(), df

def intersect_with_hierarchy(strs_bed: BedTool, region_beds: Dict[str, BedTool]) -> pl.DataFrame:
    """Perform intersections and ensure gene_id is captured"""
    all_intersections = []
    for region_name, region_bed in region_beds.items():
        if region_bed.count() == 0: continue
        
        # The BedTool command with `wb=True` retrieves the columns from file B (the region's BED file)
        # The `gene_id` is in the 4th column of file B (index 7 in the join)
        intersected = strs_bed.intersect(region_bed.sort(), wa=True, wb=True)
        
        for interval in intersected:
            # interval.fields[3] is the 'key' of the variant
            # interval.fields[7] is the 'name' of the region's BED file (which we defined as gene_id)
            g_id = interval.fields[7] if len(interval.fields) >= 8 else "."
            
            all_intersections.append({
                "key": interval.fields[3],
                "region": region_name,
                "gene_id": g_id if g_id != "." else None,
                "priority": Config.PRIORITY_MAP.get(region_name, 99)
            })
            
    if not all_intersections: return pl.DataFrame()
    
    return (pl.DataFrame(all_intersections)
            .sort(["key", "priority"])
            .unique(subset=["key"], keep="first"))

def annotate_others(df_strs: pl.DataFrame, df_annotated: pl.DataFrame) -> pl.DataFrame:
    """Mark unclassified variants as others"""
    if df_annotated.is_empty():
        return pl.DataFrame({"key": df_strs["key"], "region": "others", "gene_id": ".", "priority": Config.PRIORITY_MAP["others"]})
    
    unannotated = set(df_strs["key"]) - set(df_annotated["key"])
    if not unannotated: return df_annotated
    
    df_others = pl.DataFrame({
        "key": list(unannotated), "region": ["others"] * len(unannotated),
        "gene_id": ["."] * len(unannotated), "priority": [Config.PRIORITY_MAP["others"]] * len(unannotated)
    })
    return pl.concat([df_annotated, df_others], how="vertical")

def add_gene_information(df_annotations: pl.DataFrame, df_gtf: pl.DataFrame) -> pl.DataFrame:
    """Add complete gene metadata with strict ID mapping"""
    
    # Create a gene dictionary from the GTF (only ‘gene’ features)
    gene_index = df_gtf.filter(pl.col("feature").str.to_lowercase() == "gene").select([
        "gene_id", "gene_name", "gene_biotype"
    ]).unique(subset=["gene_id"])

    # Joining the annotations with the gene index
    df_enriched = df_annotations.join(gene_index, on="gene_id", how="left")

    # Handling of Intergenic (where gene_id must be null/empty)
    return df_enriched.with_columns([
        pl.when(pl.col("region") == "intergenic")
        .then(pl.lit("."))
        .otherwise(pl.col("gene_id").fill_null("."))
        .alias("gene_id"),
        
        pl.when(pl.col("region") == "intergenic")
        .then(pl.lit("."))
        .otherwise(pl.col("gene_name").fill_null("."))
        .alias("gene_name"),

        pl.when(pl.col("region") == "intergenic")
        .then(pl.lit("intergenic"))
        .otherwise(pl.col("gene_biotype").fill_null("unknown"))
        .alias("gene_biotype")
    ])

def calculate_statistics(df_annotations: pl.DataFrame, df_strs: pl.DataFrame) -> Dict:
    """Calculate summary statistics"""
    return {
        "summary": {"total_variants": df_strs.shape[0], "total_annotated": df_annotations.shape[0]},
        "distribution": {row["region"]: row["count"] for row in df_annotations["region"].value_counts().to_dicts()}
    }

def save_results(df_final: pl.DataFrame, stats: Dict):
    """Save output files"""
    df_final.write_csv(Path(Config.OUTDIR) / "STRs_annotated_region.tsv", separator="\t")
    stats_rows = [{"category": "summary", "metric": k, "value": v} for k, v in stats["summary"].items()]
    for region, count in stats["distribution"].items():
        stats_rows.append({"category": "distribution", "metric": region, "value": count})
    pl.DataFrame(stats_rows).write_csv(Path(Config.OUTDIR) / "code_regions_statistics.tsv", separator="\t")

# ==============================
# MAIN PIPELINE
# ==============================

def main():
    """Main function"""
    start_time = datetime.now()
    print("=" * 60)
    print("STR ANNOTATION PIPELINE (Coding vs Non-coding)")
    print("=" * 60)

    try:
        # 0. Validate input files
        print("\n[Step 0] Validating input files...")
        validate_input_files()

        # 1. Read and process GTF
        print("\n[Step 1] Reading GTF...")
        df_gtf_all = read_gtf_features(Config.GTF_PATH)

        # Filter valid chromosomes
        valid_chroms = get_valid_chromosomes(Config.GENOME_FILE)
        df_gtf_all = filter_valid_chromosomes(df_gtf_all, valid_chroms)

        if df_gtf_all.shape[0] == 0:
            print("No valid chromosomes found!")
            return

        # 2. Split Coding vs Non-Coding
        print("\n[Step 2] Splitting Coding vs Non-Coding genes...")
        df_coding = df_gtf_all.filter(pl.col("gene_biotype") == "protein_coding")
        df_non_coding = df_gtf_all.filter(pl.col("gene_biotype") != "protein_coding")

        print(f"Features: {df_coding.shape[0]} coding, {df_non_coding.shape[0]} non-coding")

        # 3. Create BED files for Hierarchical regions
        print("\n[Step 3] Generating genomic region files...")
        
        # Coding Detailed Features
        cds_bed = save_bedtool(df_coding.filter(pl.col("feature").str.to_lowercase() == "cds"), "CDS.bed")
        utr5_bed = save_bedtool(df_coding.filter(pl.col("feature").str.to_lowercase() == "five_prime_utr"), "five_prime_utr.bed")
        utr3_bed = save_bedtool(df_coding.filter(pl.col("feature").str.to_lowercase() == "three_prime_utr"), "three_prime_utr.bed")
        coding_genes_bed = save_bedtool(df_coding.filter(pl.col("feature").str.to_lowercase() == "gene"), "coding_genes.bed")
        coding_exons_bed = save_bedtool(df_coding.filter(pl.col("feature").str.to_lowercase() == "exon"), "exons.bed")
        
        promoters_bed = create_promoters(coding_genes_bed)
        introns_bed = create_introns(coding_genes_bed, coding_exons_bed)

        # Non-coding genes grouped as 'others'
        others_bed = save_bedtool(df_non_coding.filter(pl.col("feature").str.to_lowercase() == "gene"), "non_coding_genes.bed")

        # Intergenic (Complement of ALL genes)
        all_genes_merged = save_bedtool(df_gtf_all.filter(pl.col("feature").str.to_lowercase() == "gene"), "all_genes.bed").sort().merge()
        intergenic_bed = all_genes_merged.complement(g=Config.GENOME_FILE)
        intergenic_bed.saveas(str(Path(Config.OUTDIR) / "intergenic.bed"))

        # 4. Process TSV
        print("\n[Step 4] Processing TSV...")
        strs_bed, df_strs = process_tsv_to_bed(Config.TSV_PATH)

        # 5. Prepare regions for intersection
        print("\n[Step 5] Preparing regions for intersection...")
        target_regions = {
            "CDS": cds_bed,
            "five_prime_utr": utr5_bed,
            "three_prime_utr": utr3_bed,
            "non_coding_exons": coding_exons_bed,   
            "promoter": promoters_bed,
            "intron": introns_bed,
            "others": others_bed,
            "intergenic": intergenic_bed
        }

        # 6. Perform intersections
        print("\n[Step 6] Performing intersections...")
        df_annotated = intersect_with_hierarchy(strs_bed, target_regions)

        # 7. Final enrichment and saving
        df_all = annotate_others(df_strs, df_annotated)
        df_enriched = add_gene_information(df_all, df_gtf_all)

        # Keep columns from original STR analysis
        cols_to_keep = ["key", "STRs_ID", "sample_id", "chrom", "start", "end", "repeat_unit", "allele1_est", "allele2_est", "region", "gene_id", "gene_name", "depth"]
        
        df_strs_clean = df_strs.select([c for c in cols_to_keep if c in df_strs.columns])

        df_final = df_enriched.join(df_strs_clean, on="key", how="left").drop("key")
        
        # Order columns
        if "STRs_ID" in df_final.columns:
            cols = ["STRs_ID"] + [col for col in df_final.columns if col != "STRs_ID"]
            df_final = df_final.select(cols)

        # 10. Stats and Save
        stats = calculate_statistics(df_all, df_strs)
        save_results(df_final, stats)

        print(f"\nPIPELINE COMPLETED SUCCESSFULLY in {datetime.now() - start_time}")
        clean_temp_files()

    except Exception as e:
        print(f"\nFatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()