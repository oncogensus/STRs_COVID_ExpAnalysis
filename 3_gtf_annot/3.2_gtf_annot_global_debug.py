import sys
import polars as pl
from pathlib import Path
from pybedtools import BedTool


def add_chr(chrom: str) -> str:
    """Ensures the 'chr' prefix in chromosomes."""
    s_chrom = str(chrom)
    return s_chrom if s_chrom.startswith("chr") else f"chr{s_chrom}"


def run_debug():

    # =========================================================
    # Input and output file paths & Configuration
    # =========================================================
    
    # Base directory for samples
    samples_dir = Path("../samples")
    samples_dir.mkdir(parents=True, exist_ok=True)

    # Main annotated STR dataset
    results_path = samples_dir / "STRs_annotated_region.tsv"

    # Output file that will store all variants classified as "others"
    output_others_path = samples_dir / "others_regions_statistics.csv"
    
    # Audit configuration
    gtf_path = Path("Homo_sapiens.GRCh38.98.gtf")
    
    # =========================================================
    # Check if input file exists before continuing
    # =========================================================
    if not results_path.exists():
        print(f"Error: File {results_path} not found.")
        return

    print("=== STR ANNOTATION DEBUG LOG ===")

    try:
        # =========================================================
        # Read TSV file using Polars
        # =========================================================
        df = pl.read_csv(
            results_path,
            separator="\t",
            infer_schema_length=10000,
            schema_overrides={
                "chrom": pl.String,
                "gene_id": pl.String,
                "STRs_ID": pl.String,
                "gene_biotype": pl.String
            }
        )

    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # =========================================================
    # [1] Global distribution of variants by genomic region
    # =========================================================
    print("\n[1] Distribution of Variants by Region (coding region):")

    if "region" in df.columns:
        dist = (
            df["region"]
            .value_counts()
            .sort("count", descending=True)
        )
        print(dist)
    else:
        print("Column 'region' not found!")

    # =========================================================
    # [2] Detailed analysis of variants classified as "others"
    # =========================================================
    print("\n[2] Detailed Analysis: 'others' regions")

    df_others = df.filter(pl.col("region") == "others")

    if df_others.is_empty():
        print(" - No variants found in 'others' region.")
        return # Encerrar pois não há "others" para auditar

    total_others = df_others.shape[0]
    print(f" - Total 'others' variants found: {total_others}")

    if "gene_biotype" in df_others.columns:
        biotype_dist = (
            df_others["gene_biotype"]
            .value_counts()
            .sort("count", descending=True)
        )
        print("\nBiotype distribution within 'others':")
        print(biotype_dist)

    try:
        df_others.write_csv(output_others_path)
        print(f"\n -> SUCCESS: 'others' location report saved to: {output_others_path}")
    except Exception as e:
        print(f" -> ERROR saving others_loc.csv: {e}")

    # =========================================================
    # [3] Integrity check
    # =========================================================
    print("\n[3] Integrity Check (others vs intergenic):")

    if all(col in df.columns for col in ["region", "gene_id"]):
        others_with_gene = df.filter(
            (pl.col("region") == "others") &
            (pl.col("gene_id") != ".")
        ).shape[0]

        intergenic_with_gene = df.filter(
            (pl.col("region") == "intergenic") &
            (pl.col("gene_id") != ".")
        ).shape[0]

        print(f" - 'others' variants with Gene ID: {others_with_gene}/{df_others.shape[0]}")
        print(f" - 'intergenic' variants with Gene ID: {intergenic_with_gene} (Expected: 0)")

    # =========================================================
    # [4] Display a small sample of "others" variants
    # =========================================================
    print("\n[4] 'Others' Variants Sample (Top 5):")

    cols_view = ["STRs_ID", "chrom", "start", "end", "gene_name", "gene_biotype"]
    available_cols = [c for c in cols_view if c in df_others.columns]

    print(df_others.select(available_cols).head(5))

    # =========================================================
    # [5] Protein Coding Audit for 'others'
    # =========================================================
    print("\n[5] STR 'Others' - Protein Coding Audit:")
    
    if not gtf_path.exists():
        print(f" - SKIPPING AUDIT: GTF file {gtf_path} not found.")
    else:
        print(" - Formatting coordinates for BedTool...")
        
        # Format df_others for BED intersection
        bed_df_others = df_others.with_columns([
            pl.col("chrom").map_elements(add_chr, return_dtype=pl.String),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64)
        ])
        
        # Salvando na pasta ../samples/
        bed_path = samples_dir / "temp_others.bed"
        bed_df_others.select([
            "chrom", 
            (pl.col("start") - 1).alias("start"), 
            "end", 
            "STRs_ID"
        ]).write_csv(bed_path, separator="\t", include_header=False)
        
        others_bed = BedTool(str(bed_path)).sort()

        print(f" - Processing GTF for protein_coding features...")
        gtf_columns = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        
        # Read GTF
        df_gtf = pl.read_csv(
            gtf_path, 
            separator="\t", 
            has_header=False, 
            comment_prefix="#", 
            new_columns=gtf_columns,
            infer_schema_length=0 
        )
        
        # Filter protein_coding and process coordinates
        df_coding = df_gtf.filter(
            pl.col("attribute").str.contains('gene_biotype "protein_coding"')
        ).with_columns([
            pl.col("chrom").map_elements(add_chr, return_dtype=pl.String),
            pl.col("attribute").str.extract(r'gene_name "([^"]+)"', 1).alias("gene_name"),
            pl.col("start").cast(pl.Int64),
            pl.col("end").cast(pl.Int64)
        ])

        # Salvando na pasta ../samples/
        gtf_bed_path = samples_dir / "temp_coding.bed"
        df_coding.select([
            "chrom", 
            (pl.col("start") - 1).alias("start"), 
            "end", 
            (pl.col("feature") + "|" + pl.col("gene_name").fill_null("Unknown")).alias("metadata")
        ]).write_csv(gtf_bed_path, separator="\t", include_header=False)
        
        coding_bed = BedTool(str(gtf_bed_path)).sort()

        # Intersection
        print(" - Performing intersection to identify overlap origins...")
        intersected = others_bed.intersect(coding_bed, wa=True, wb=True)
        
        records = []
        for interval in intersected:
            metadata = interval.fields[7].split("|")
            records.append({
                "STRs_ID": interval.fields[3],
                "overlap_feature": metadata[0],
                "gene_name": metadata[1]
            })

        if not records:
            print(" - No protein_coding overlaps found within the 'others' category.")
        else:
            df_audit = pl.DataFrame(records)
            print("\n === AUDIT RESULTS: Protein Coding Overlaps in 'Others' ===")
            
            unique_strs = df_audit["STRs_ID"].n_unique()
            print(f" Total unique 'others' variants overlapping protein_coding: {unique_strs}")

            feature_dist = df_audit.group_by("overlap_feature").len().sort("len", descending=True)
            print("\n Distribution by Feature Type:")
            print(feature_dist)

        # Cleanup temporary files (remove da pasta ../samples/)
        if bed_path.exists(): bed_path.unlink()
        if gtf_bed_path.exists(): gtf_bed_path.unlink()

    # =========================================================
    # End of debug routine
    # =========================================================
    print("\n" + "=" * 30)
    print("DEBUG COMPLETED")


if __name__ == "__main__":
    run_debug()