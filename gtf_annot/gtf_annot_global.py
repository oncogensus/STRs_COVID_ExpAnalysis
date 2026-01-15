"""
STR Annotation Pipeline - Python Version
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
    # Input files
    TSV_PATH = "../samples/gp_global/merged/str_consolidated.tsv"
    GTF_PATH = "Homo_sapiens.GRCh38.98.gtf"
    GENOME_FILE = "genome.txt"

    # Output directory
    OUTDIR = "../samples/gp_global/merged/"

    # Parameters
    PROMOTER_DIST = 3000
    
    # Annotation priorities
    PRIORITY_MAP = {
        "CDS": 1,
        "five_prime_utr": 2,
        "three_prime_utr": 3,
        "promoter": 4,
        "intron": 5,
        "intergenic": 6,
        "others": 7
    }

    # Columns
    GENE_ID_COLUMN = "gene_id"
    KEEP_TEMP_FILES = False

# Create output directory
Path(Config.OUTDIR).mkdir(parents=True, exist_ok=True)

# ==============================
# HELPER FUNCTIONS
# ==============================

def add_chr(chrom: str) -> str:
    """Add 'chr' to chromosome if not present"""
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"

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
            "genes.bed", "exons.bed", "CDS.bed",
            "five_prime_utr.bed", "three_prime_utr.bed",
            "promoters.bed", "introns.bed", "intergenic.bed"
        ]
        for file in temp_files:
            path = Path(Config.OUTDIR) / file
            if path.exists():
                try:
                    path.unlink()
                except:
                    pass

def read_gtf_features(gtf_path: str) -> pl.DataFrame:
    """Read GTF file and extract attributes"""
    print(f"Reading GTF: {gtf_path}")

    # Correct schema with 9 columns
    schema = {
        "chrom": pl.Utf8,
        "source": pl.Utf8,
        "feature": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "frame": pl.Utf8,
        "attribute": pl.Utf8
    }

    # Read GTF
    df = pl.read_csv(
        gtf_path,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        null_values=".",
        new_columns=list(schema.keys()),
        schema_overrides=schema
    )

    # Extract gene_id and gene_name
    df = df.with_columns([
        pl.col("attribute").str.extract(r'gene_id\s+"([^"]+)"', 1).alias("gene_id"),
        pl.col("attribute").str.extract(r'gene_name\s+"([^"]+)"', 1).alias("gene_name")
    ])

    # Use gene_id as fallback for gene_name
    df = df.with_columns(
        pl.when(pl.col("gene_name").is_null())
        .then(pl.col("gene_id"))
        .otherwise(pl.col("gene_name"))
        .alias("gene_name")
    )

    # Add 'chr' and select columns
    df = df.with_columns(
        pl.col("chrom").map_elements(add_chr, return_dtype=pl.Utf8).alias("chrom")
    ).select([
        "chrom", "start", "end", "feature",
        "gene_id", "gene_name"
    ])

    print(f"GTF read: {df.shape[0]} rows")
    return df

def get_valid_chromosomes(genome_file: str) -> Set[str]:
    """Read genome file and return valid chromosomes"""
    valid_chroms = set()
    with open(genome_file, 'r') as f:
        for line in f:
            if line.strip():
                chrom = line.split('\t')[0].strip()
                valid_chroms.add(chrom)
    print(f"Valid chromosomes: {len(valid_chroms)}")
    return valid_chroms

def filter_valid_chromosomes(df: pl.DataFrame, valid_chroms: Set[str]) -> pl.DataFrame:
    """Filter DataFrame to keep only valid chromosomes"""
    original_count = df.shape[0]
    df_filtered = df.filter(pl.col("chrom").is_in(list(valid_chroms)))

    if df_filtered.shape[0] < original_count:
        removed = set(df.filter(~pl.col("chrom").is_in(list(valid_chroms)))["chrom"].unique().to_list())
        print(f"Removed {original_count - df_filtered.shape[0]} records from {len(removed)} invalid chromosomes")

    return df_filtered

def save_bedtool(df: pl.DataFrame, filename: str, name_col: str = None) -> BedTool:
    """Save DataFrame as BED file and return BedTool - 4 COLUMNS ONLY"""
    path = Path(Config.OUTDIR) / filename

    # Use default column if not specified
    if name_col is None:
        name_col = Config.GENE_ID_COLUMN

    # Check if column exists
    if name_col not in df.columns:
        print(f"Warning: Column {name_col} not found. Using '.'")
        df = df.with_columns(pl.lit(".").alias("name"))
        name_col = "name"

    # Convert coordinates 1-based → 0-based
    bed_df = df.with_columns([
        (pl.col("start") - 1).alias("start"),
        pl.col("end").alias("end")
    ]).select([
        "chrom", "start", "end",
        pl.col(name_col).alias("name")
    ])

    # Save
    bed_df.write_csv(path, separator="\t", include_header=False)
    print(f"Saved: {path} ({bed_df.shape[0]} records)")

    return BedTool(str(path))

def create_gene_regions(df_gtf: pl.DataFrame) -> Dict[str, BedTool]:
    """Create basic genomic regions"""
    print("Creating gene regions...")

    features = {}

    # Genes
    genes = df_gtf.filter(pl.col("feature").str.to_lowercase() == "gene")
    features["gene"] = save_bedtool(genes, "genes.bed")

    # Exons
    exons = df_gtf.filter(pl.col("feature").str.to_lowercase() == "exon")
    features["exon"] = save_bedtool(exons, "exons.bed")

    # CDS
    cds = df_gtf.filter(pl.col("feature").str.to_lowercase() == "cds")
    features["CDS"] = save_bedtool(cds, "CDS.bed")

    # UTRs
    utr5 = df_gtf.filter(pl.col("feature").str.to_lowercase() == "five_prime_utr")
    features["five_prime_utr"] = save_bedtool(utr5, "five_prime_utr.bed")

    utr3 = df_gtf.filter(pl.col("feature").str.to_lowercase() == "three_prime_utr")
    features["three_prime_utr"] = save_bedtool(utr3, "three_prime_utr.bed")

    return features

def create_promoters(genes_bed: BedTool) -> BedTool:
    """Create promoter regions upstream of TSS"""
    print("Creating promoters...")

    if genes_bed.count() == 0:
        print("Warning: No genes for promoter creation")
        return BedTool("", from_string=True)

    # Create TSS (Transcription Start Site)
    tss_records = []
    for gene in genes_bed:
        tss_start = gene.start
        tss_end = tss_start + 1

        tss_records.append([
            gene.chrom, tss_start, tss_end,
            gene.name
        ])

    if not tss_records:
        return BedTool("", from_string=True)

    tss_bed = BedTool("\n".join("\t".join(map(str, r)) for r in tss_records), from_string=True)

    # Expand to promoter region
    promoters = tss_bed.slop(
        g=Config.GENOME_FILE,
        l=Config.PROMOTER_DIST,
        r=0,
        s=False
    )

    # Save
    promoters.saveas(str(Path(Config.OUTDIR) / "promoters.bed"))
    print(f"Promoters created: {len(promoters)}")
    return promoters

def create_introns(genes_bed: BedTool, exons_bed: BedTool) -> BedTool:
    """Create introns by subtracting exons from genes"""
    print("Creating introns...")
    
    if genes_bed.count() == 0:
        print("Warning: No genes for intron creation")
        return BedTool("", from_string=True)

    # Group exons by gene
    exons_by_gene = defaultdict(list)
    for exon in exons_bed:
        exons_by_gene[exon.name].append(exon)

    # Process each gene
    intron_records = []
    for gene in genes_bed:
        gene_id = gene.name

        if gene_id in exons_by_gene:
            gene_exons = BedTool(exons_by_gene[gene_id]).sort()
            gene_bed = BedTool(f"{gene.chrom}\t{gene.start}\t{gene.end}\t{gene.name}", from_string=True)
            
            # Subtract exons to get introns
            introns = gene_bed.subtract(gene_exons)

            for intron in introns:
                intron_records.append([
                    intron.chrom, intron.start, intron.end,
                    gene_id
                ])
        else:
            # Gene without annotated exons
            intron_records.append([
                gene.chrom, gene.start, gene.end,
                gene_id
            ])

    # Create BedTool
    if intron_records:
        introns_bed = BedTool("\n".join("\t".join(map(str, r)) for r in intron_records), from_string=True).sort()
        print(f"Introns created: {len(intron_records)}")
    else:
        introns_bed = BedTool("", from_string=True)
        print("Warning: No introns created")

    # Save
    introns_bed.saveas(str(Path(Config.OUTDIR) / "introns.bed"))
    return introns_bed

def create_intergenic_regions(genes_bed: BedTool) -> BedTool:
    """Create intergenic regions"""
    print("Creating intergenic regions...")
    
    if genes_bed.count() == 0:
        print("Warning: No genes for intergenic region creation")
        return BedTool("", from_string=True)

    # Merge genes
    genes_merged = genes_bed.sort().merge()

    # Complement to get intergenic regions
    intergenic = genes_merged.complement(g=Config.GENOME_FILE)

    # Add metadata
    intergenic_records = []
    for region in intergenic:
        intergenic_records.append([
            region.chrom, region.start, region.end,
            "intergenic"
        ])

    intergenic_bed = BedTool("\n".join("\t".join(map(str, r)) for r in intergenic_records), from_string=True)

    # Save
    intergenic_bed.saveas(str(Path(Config.OUTDIR) / "intergenic.bed"))
    print(f"Intergenic regions: {len(intergenic)}")
    return intergenic_bed

def process_tsv_to_bed(tsv_path: str) -> Tuple[BedTool, pl.DataFrame]:
    """Process consolidated TSV file"""
    print(f"Processing TSV: {tsv_path}")
    
    try:
        # Read TSV with tab separator
        df_tsv = pl.read_csv(
            tsv_path,
            separator="\t",
            has_header=True,
            null_values=["", "NA", "null", "NaN", ".", "-", "None", "NULL"],
            infer_schema_length=10000,
            try_parse_dates=False
        )
        
        print(f"TSV rows: {df_tsv.shape[0]}")
        print(f"Columns: {df_tsv.columns}")
        
    except Exception as e:
        print(f"Error reading TSV: {e}")
        import traceback
        traceback.print_exc()
        raise

    # Check required columns
    required_cols = ["chrom", "start", "end", "repeat_unit"]
    for col in required_cols:
        if col not in df_tsv.columns:
            raise ValueError(f"Column {col} not found. Available: {df_tsv.columns}")

    # Filter null values
    df_filtered = df_tsv.filter(
        pl.col("chrom").is_not_null() &
        pl.col("start").is_not_null() &
        pl.col("end").is_not_null() &
        pl.col("repeat_unit").is_not_null()
    )
    
    if df_filtered.shape[0] == 0:
        print("Error: No valid rows found!")
        print(f"Total rows: {df_tsv.shape[0]}")
        print("Null counts:")
        for col in required_cols:
            null_count = df_tsv[col].null_count()
            print(f"  {col}: {null_count} nulls")
        return BedTool("", from_string=True), pl.DataFrame()
    
    # Add unique ID and process coordinates - SEM DUPLICATAS
    df = df_filtered.with_row_index("row_id").with_columns([
        # Add 'chr' if necessary
        pl.when(pl.col("chrom").str.starts_with("chr"))
        .then(pl.col("chrom"))
        .otherwise("chr" + pl.col("chrom"))
        .alias("var_chrom"),
        
        # Convert start to 0-based for BED
        (pl.col("start") - 1).alias("var_start"),
        
        # End remains the same (BED is exclusive)
        pl.col("end").alias("var_end"),
        
        # Create unique key (temporário para junções)
        pl.col("row_id").cast(pl.Utf8).alias("key"),
        
        # Calculate length
        (pl.col("end") - pl.col("start") + 1).alias("region_length"),
        
        # Keep sample if exists
        pl.col("sample").alias("sample_id") if "sample" in df_filtered.columns
        else pl.lit("unknown").alias("sample_id")
    ])
    
    # Calculate repeat count (assuming region_length is multiple of repeat_unit length)
    df = df.with_columns([
        (pl.col("region_length") / pl.col("repeat_unit").str.len_bytes()).alias("repeat_count")
    ])
    
    # Create STRs_ID column
    df = df.with_columns([
        (
            pl.col("var_chrom") + ":" + 
            pl.col("start").cast(pl.Utf8) + ":" + 
            pl.col("repeat_unit") + ":" + 
            pl.col("repeat_count").cast(pl.Int64).cast(pl.Utf8)
        ).alias("STRs_ID")
    ])
    
    # Create BED file
    if df.shape[0] > 0:
        bed_path = Path(Config.OUTDIR) / "str_regions_temp.bed"
        
        # Create BED DataFrame - only 4 columns needed for bedtools
        bed_df = df.select([
            "var_chrom", "var_start", "var_end", "key"
        ])
        
        bed_df.write_csv(
            bed_path,
            separator="\t",
            include_header=False
        )
        
        bed = BedTool(str(bed_path)).sort()
        
        if not Config.KEEP_TEMP_FILES:
            bed_path.unlink()
        
        print(f"\nCreated BED with {bed.count()} STR regions")
        print(f"  Format: 4 columns (chrom, start, end, name)")
        print(f"  Sample chromosomes: {df['var_chrom'].head(3).to_list()}")
    else:
        bed = BedTool("", from_string=True)
    
    return bed, df

def intersect_with_hierarchy(strs_bed: BedTool, region_beds: Dict[str, BedTool]) -> pl.DataFrame:
    """Perform intersections with priority hierarchy"""
    print("Performing intersections...")
    print(f"Total STR regions: {strs_bed.count()}")
    
    all_intersections = []

    for region_name, region_bed in region_beds.items():
        if region_bed.count() == 0:
            print(f"  Skipping {region_name}: empty region")
            continue

        print(f"  Intersecting with {region_name}...")
        print(f"    STR regions: {strs_bed.count()}")
        print(f"    {region_name} regions: {region_bed.count()}")

        try:
            intersected = strs_bed.intersect(region_bed.sort(), wa=True, wb=True)
        except Exception as e:
            print(f"    Error in intersection: {e}")
            continue

        # Process results
        intersection_count = 0
        for interval in intersected:
            str_key = interval.fields[3]
            
            gene_id = None
            if len(interval.fields) >= 8:
                gene_id = interval.fields[7] if interval.fields[7] != "." else None
            
            all_intersections.append({
                "key": str_key,
                "region": region_name,
                "gene_id": gene_id,
                "priority": Config.PRIORITY_MAP.get(region_name, 999)
            })
            intersection_count += 1
        
        print(f"    Found {intersection_count} intersections")

    # Return DataFrame with proper schema
    if not all_intersections:
        return pl.DataFrame({
            "key": pl.Series([], dtype=pl.Utf8),
            "region": pl.Series([], dtype=pl.Utf8),
            "gene_id": pl.Series([], dtype=pl.Utf8),
            "priority": pl.Series([], dtype=pl.Int64)
        })

    df = pl.DataFrame(all_intersections)
    
    # Apply hierarchy
    df = df.sort(["key", "priority"]).unique(subset=["key"], keep="first")

    print(f"\nIntersections completed: {df.shape[0]} annotated variants")
    return df

def add_gene_information(df_annotations: pl.DataFrame, df_gtf: pl.DataFrame) -> pl.DataFrame:
    """Add complete gene information"""
    print("Adding gene information...")

    # Create gene index
    gene_index = df_gtf.filter(
        pl.col("feature").str.to_lowercase() == "gene"
    ).select([
        "gene_id", "gene_name", "chrom", "start", "end"
    ]).unique(subset=["gene_id"], keep="first").rename({
        "chrom": "gene_chrom",
        "start": "gene_start",
        "end": "gene_end"
    })

    # Combine with annotations
    df_enriched = df_annotations.join(gene_index, on="gene_id", how="left")
    
    # Fill null values
    df_enriched = df_enriched.with_columns([
        pl.col("gene_id").fill_null("."),
        pl.col("gene_name").fill_null("."),
        pl.col("gene_chrom").fill_null("."),
        pl.col("gene_start").fill_null(-1),
        pl.col("gene_end").fill_null(-1)
    ])

    # Create combined annotation
    df_enriched = df_enriched.with_columns(
        pl.when(pl.col("gene_name") == ".")
        .then(pl.col("region"))
        .otherwise(pl.col("region") + ":" + pl.col("gene_name"))
        .alias("annotation")
    )

    return df_enriched

def annotate_others(df_strs: pl.DataFrame, df_annotated: pl.DataFrame) -> pl.DataFrame:
    """Identify unclassified variants"""
    all_keys = set(df_strs["key"].to_list())
    annotated_keys = set(df_annotated["key"].to_list())
    unannotated_keys = all_keys - annotated_keys

    print(f"Total variants: {len(all_keys)}")
    print(f"Annotated: {len(annotated_keys)}")
    print(f"Unannotated: {len(unannotated_keys)}")

    if unannotated_keys:
        df_others = pl.DataFrame({
            "key": list(unannotated_keys),
            "region": ["others"] * len(unannotated_keys),
            "gene_id": ["."] * len(unannotated_keys),
            "priority": [Config.PRIORITY_MAP["others"]] * len(unannotated_keys)
        })
        
        if df_annotated.is_empty():
            return df_others
        
        df_annotated = df_annotated.with_columns([
            pl.col("key").cast(pl.Utf8),
            pl.col("region").cast(pl.Utf8),
            pl.col("gene_id").cast(pl.Utf8),
            pl.col("priority").cast(pl.Int64)
        ])
        
        return pl.concat([df_annotated, df_others], how="vertical")

    return df_annotated

def calculate_statistics(df_annotations: pl.DataFrame, df_strs: pl.DataFrame) -> Dict:
    """Calculate annotation statistics"""
    print("Calculating statistics...")

    stats = {
        "summary": {
            "total_variants": df_strs.shape[0],
            "total_annotated": df_annotations.shape[0]
        }
    }

    if "region" in df_annotations.columns:
        region_counts = df_annotations["region"].value_counts().to_dicts()
        stats["distribution"] = {row["region"]: row["count"] for row in region_counts}

    return stats

def save_results(df_final: pl.DataFrame, stats: Dict):
    """Save results"""
    print("Saving results...")

    # Main file
    output_file = Path(Config.OUTDIR) / "STRs_annotated_complete.tsv"
    df_final.write_csv(output_file, separator="\t")
    print(f"Main annotations: {output_file}")

    # Statistics
    stats_file = Path(Config.OUTDIR) / "global_annotation_statistics.tsv"
    stats_rows = [
        {"category": "summary", "metric": "total_variants", "value": stats["summary"]["total_variants"]},
        {"category": "summary", "metric": "total_annotated", "value": stats["summary"]["total_annotated"]}
    ]

    if "distribution" in stats:
        for region, count in stats["distribution"].items():
            stats_rows.append({
                "category": "distribution",
                "metric": region,
                "value": count
            })

    pl.DataFrame(stats_rows).write_csv(stats_file, separator="\t")
    print(f"Statistics: {stats_file}")

def main():
    """Main function"""
    start_time = datetime.now()
    print("=" * 60)
    print("STR ANNOTATION PIPELINE")
    print("=" * 60)

    try:
        # 0. Validate input files
        print("\n[Step 0] Validating input files...")
        validate_input_files()

        # 1. Read and process GTF
        print("\n[Step 1] Reading GTF...")
        df_gtf = read_gtf_features(Config.GTF_PATH)

        # Filter valid chromosomes
        valid_chroms = get_valid_chromosomes(Config.GENOME_FILE)
        df_gtf = filter_valid_chromosomes(df_gtf, valid_chroms)

        if df_gtf.shape[0] == 0:
            print("No valid chromosomes found!")
            return

        # 2. Create gene regions
        print("\n[Step 2] Creating gene regions...")
        basic_regions = create_gene_regions(df_gtf)

        # 3. Create special regions
        print("\n[Step 3] Creating special regions...")
        promoters = create_promoters(basic_regions["gene"])
        introns = create_introns(basic_regions["gene"], basic_regions["exon"])
        intergenic = create_intergenic_regions(basic_regions["gene"])

        # 4. Process TSV
        print("\n[Step 4] Processing TSV...")
        strs_bed, df_strs = process_tsv_to_bed(Config.TSV_PATH)

        if strs_bed.count() == 0:
            print("No variants found!")
            return

        # 5. Prepare regions for intersection
        print("\n[Step 5] Preparing regions for intersection...")
        all_regions = {
            "CDS": basic_regions["CDS"],
            "five_prime_utr": basic_regions["five_prime_utr"],
            "three_prime_utr": basic_regions["three_prime_utr"],
            "promoter": promoters,
            "intron": introns,
            "intergenic": intergenic
        }

        # 6. Perform intersections
        print("\n[Step 6] Performing intersections...")
        df_annotated = intersect_with_hierarchy(strs_bed, all_regions)

        # 7. Add "others"
        print("\n[Step 7] Identifying unclassified variants...")
        df_all = annotate_others(df_strs, df_annotated)

        # 8. Add gene information
        print("\n[Step 8] Adding gene information...")
        df_enriched = add_gene_information(df_all, df_gtf)

        # Keep only necessary columns from df_strs
        cols_to_keep = [
            "key", "STRs_ID", "sample_id", "chrom", "start", "end", "repeat_unit",
            "allele1_est", "allele2_est", "depth"
        ]

        df_strs_clean = df_strs.select(cols_to_keep)

        # 9. Combine with original information
        print("\n[Step 9] Combining with original data...")
        df_final = df_enriched.join(
            df_strs_clean,
            on="key",
            how="left"
        )
        
        # Remove 'key' column and reorder columns
        df_final = df_final.drop("key")
        
        # Reorder columns to put STRs_ID first
        cols = df_final.columns
        if "STRs_ID" in cols:
            # Move STRs_ID to the beginning
            cols = ["STRs_ID"] + [col for col in cols if col != "STRs_ID"]
            df_final = df_final.select(cols)
        
        # 10. Calculate statistics and save
        print("\n[Step 10] Calculating statistics...")
        stats = calculate_statistics(df_all, df_strs)

        print("\n[Step 11] Saving results...")
        save_results(df_final, stats)

        # Final summary
        elapsed = datetime.now() - start_time
        print("\n" + "=" * 60)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print(f"Total time: {elapsed}")
        print(f"Variants processed: {stats['summary']['total_variants']}")
        print(f"Variants annotated: {stats['summary']['total_annotated']}")

        if "distribution" in stats:
            print("\nDistribution by region:")
            for region, count in stats["distribution"].items():
                percentage = (count / stats['summary']['total_annotated'] * 100) if stats['summary']['total_annotated'] > 0 else 0
                print(f"  {region}: {count} ({percentage:.1f}%)")

        print("=" * 60)

        # Clean temporary files
        clean_temp_files()

    except Exception as e:
        print(f"\nFatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()