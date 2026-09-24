import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
from great_tables import GT, html

BASE = Path(__file__).resolve().parent

PT_TO_EN = {
    "Pulmão e Cérebro": "Lung and Brain",
    "Pulmão": "Lung",
    "Cérebro": "Brain",
}

GENEHANCER_IDS = {
    "GNG7": "GH19J002518 / GH19J002522 / GH19J002539",
}
GENEHANCER_SCORES = {
    "GNG7": "4 (5.640) / 5 (9.320) / 2 (2.360)",
}

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: parse_outlier
# ─────────────────────────────────────────────────────────────────────────────
def parse_outlier(html_path):
    """Return {gene: {Region, Chr, Start, Motif}} from outlier_detail_table.html."""
    soup = BeautifulSoup(Path(html_path).read_text(encoding="utf-8"), "html.parser")
    variant_map = {}
    for tr in soup.find("tbody").find_all("tr"):
        if "gt_group_heading_row" in tr.get("class", []):
            continue
        cells = tr.find_all("td")
        if len(cells) < 5:
            continue
        gene = cells[0].get_text(strip=True)
        if gene not in variant_map:
            variant_map[gene] = {
                "Region": cells[1].get_text(strip=True),
                "Chr": cells[2].get_text(strip=True),
                "Start": cells[3].get_text(strip=True),
                "Motif": cells[4].get_text(strip=True),
            }
    return variant_map

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: translate_pt
# ─────────────────────────────────────────────────────────────────────────────
def translate_pt(val):
    for pt, en in PT_TO_EN.items():
        val = val.replace(pt, en)
    return val

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: split_gene_value
# ─────────────────────────────────────────────────────────────────────────────
def split_gene_value(raw, gene_name):
    raw = str(raw).strip()
    if raw in ("nan", "—"):
        return "—"
    raw = translate_pt(raw)
    if gene_name + ":" not in raw:
        return raw
    pattern = re.compile(
        rf'{re.escape(gene_name)}:\s*([^,]*?)(?=(?:\b(?:CAMK4|DRAIC)\b:|$))',
        re.IGNORECASE,
    )
    m = pattern.search(raw)
    return m.group(1).strip() if m else raw

def format_tissue_values(raw):
    """Convert tissue values to 'Lung and Brain / Lung' format."""
    raw = str(raw).strip()
    if raw in ("nan", "—"):
        return "—"
    # Split by newlines and clean up
    lines = [l.strip() for l in raw.split("\n") if l.strip()]
    formatted = []
    for line in lines:
        # Remove tissue prefix and keep values only
        if ":" in line:
            tissue, values = line.split(":", 1)
            tissue = tissue.strip()
            values = values.strip()
            # Convert "Lung and Brain" to shorter form
            if "and" in tissue:
                tissue = "Lung/Brain"
            formatted.append(f"{tissue}: {values}")
        else:
            formatted.append(line)
    # Join with " / " if multiple lines
    return " / ".join(formatted) if formatted else "—"

def format_encode_value(raw):
    """Convert 'Max (Min)' to 'min - max' format."""
    raw = str(raw).strip()
    if raw in ("nan", "—"):
        return "—"
    import re
    match = re.match(r'([\d.]+)\s*\(([\d.]+)\)', raw)
    if match:
        max_val, min_val = match.group(1), match.group(2)
        return f"{min_val} - {max_val}"
    return raw

def format_tissue_encode(raw):
    """Format tissue-specific ENCODE values with tissue labels."""
    raw = str(raw).strip()
    if raw in ("nan", "—"):
        return "—"
    # Split by " / Lung:" to separate tissue-specific values
    if " / Lung:" in raw:
        parts = raw.split(" / Lung:")
        first_set = parts[0].strip()
        second_set = parts[1].strip() if len(parts) > 1 else ""
        # Parse first set (Lung/Brain)
        if ":" in first_set:
            _, values = first_set.split(":", 1)
            values = values.strip()
        else:
            values = first_set
        # Parse second set (Lung)
        if second_set and ":" in second_set:
            _, lung_values = second_set.split(":", 1)
            lung_values = lung_values.strip()
        else:
            lung_values = ""
        # Convert to min - max format
        first_val = format_encode_value(values.split(" / ")[0] if " / " in values else values)
        second_val = format_encode_value(values.split(" / ")[1] if " / " in values and len(values.split(" / ")) > 1 else "")
        lung_first = format_encode_value(lung_values.split(" / ")[0] if " / " in lung_values else lung_values)
        lung_second = format_encode_value(lung_values.split(" / ")[1] if " / " in lung_values and len(lung_values.split(" / ")) > 1 else "")
        return f"Lung/Brain: {first_val} / Lung: {lung_first}"
    elif " / " in raw:
        parts = raw.split(" / ")
        return format_encode_value(parts[0])
    else:
        return format_encode_value(raw)
    import re
    match = re.match(r'([\d.]+)\s*\(([\d.]+)\)', raw)
    if match:
        max_val, min_val = match.group(1), match.group(2)
        return f"{min_val} - {max_val}"
    return raw

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: build_dataframe
# ─────────────────────────────────────────────────────────────────────────────
def read_input_csv(path):
    """Read input CSV regardless of encoding (utf-8 first, fallback to cp1252/latin-1)."""
    for enc in ("utf-8", "cp1252"):
        try:
            return pd.read_csv(path, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, encoding="latin-1")

def build_dataframe(csv_path, variant_map):
    df_tp = read_input_csv(csv_path)
    records = []
    for _, row in df_tp.iterrows():
        locus = row.iloc[0]
        genes = [g.strip() for g in re.split(r"[/,]", locus)]
        coord_str, eh_rm, ccre = row.iloc[1], row.iloc[2], row.iloc[3]
        dnase_atac, histone_ctcf, jarvis = str(row.iloc[4]), str(row.iloc[5]), str(row.iloc[6])
        
        # Parse STR structure from coordinates (e.g., "chr4:... · (TAGA)14 (1.00)")
        motif_ext = "—"
        purity = "—"
        if "·" in str(coord_str):
            parts = str(coord_str).split("·")
            trexplorer = parts[1].strip() if len(parts) > 1 else "—"
            # Split "(TAGA)14 (1.00)" into motif_ext and purity
            if "(" in trexplorer:
                # Handle multiple motifs like "(AT)10 (1.00) + (AC)4 (1.00)"
                motif_parts = re.findall(r'\(([^)]+)\)\d+', trexplorer)
                purity_parts = re.findall(r'\(1\.0\d*\)', trexplorer)
                if motif_parts and purity_parts:
                    motif_ext = " / ".join([f"({m})" for m in motif_parts])
                    # Get the number of repeats
                    repeat_nums = re.findall(r'\)(\d+)', trexplorer)
                    if repeat_nums:
                        motif_ext = " / ".join([f"({motif_parts[i]}){repeat_nums[i]}" for i in range(len(motif_parts))])
                    purity = purity_parts[0]
        
        # Parse ExpHet and RM from eh_rm (e.g., "0.890 / Simple_repeat (0.0%)")
        exp_het = "—"
        rm_div = "—"
        if " / " in str(eh_rm):
            eh_parts = str(eh_rm).split(" / ", 1)
            exp_het = eh_parts[0].strip()
            rm_raw = eh_parts[1].strip() if len(eh_parts) > 1 else "—"
            # Handle cases like "0.702 · Simple_repeat (20.9%)" - extract RM part
            if "·" in rm_raw:
                rm_parts = rm_raw.split("·", 1)
                exp_het = exp_het + " / " + rm_parts[0].strip()
                rm_div = rm_parts[1].strip()
            else:
                rm_div = rm_raw
            # Format RM class name for publication (e.g., "Simple_repeat" → "Simple repeat")
            rm_div = rm_div.replace("Simple_repeat", "Simple repeat")
        
        for g in genes:
            v = variant_map.get(g, {})
            jarvis_raw = split_gene_value(jarvis, g)
            # Split DepletionRank from JARVIS if present
            if "DepletionRank" in jarvis_raw:
                # Handle "0.931 ± 0.015 (DepletionRank)" format
                depletion_match = re.search(r'([\d.]+\s*±\s*[\d.]+)\s*\(DepletionRank\)', jarvis_raw)
                if depletion_match:
                    depletion = depletion_match.group(1)
                    jarvis_clean = "—"
                else:
                    parts = jarvis_raw.split("DepletionRank")
                    jarvis_clean = parts[0].strip().rstrip(",").strip()
                    depletion = "DepletionRank" + parts[1] if len(parts) > 1 else "—"
            else:
                jarvis_clean = jarvis_raw
                depletion = "—"
            
            # Format JARVIS as Min - Max
            if jarvis_clean != "—" and "(" in jarvis_clean:
                m = re.match(r'([\d.eE+-]+)\s*\(([\d.eE+-]+)\)', jarvis_clean)
                if m:
                    jarvis_clean = f"{m.group(2)} - {m.group(1)}"
                else:
                    jarvis_clean = jarvis_clean.replace("(", "- ").replace(")", "")
            
            # Parse histone/CTCF values - handle tissue-specific (multi-line) or single-line
            histone_raw = str(split_gene_value(histone_ctcf, g)).strip()
            h3k4me3 = "—"
            h3k27ac = "—"
            ctcf = "—"
            if histone_raw not in ("nan", "—", ""):
                lines = [l.strip() for l in histone_raw.split("\n") if l.strip()]
                if len(lines) >= 2:
                    t1 = translate_pt(lines[0].split(":", 1)[0].strip())
                    t2 = translate_pt(lines[1].split(":", 1)[0].strip())
                    _, vals1 = lines[0].split(":", 1)
                    _, vals2 = lines[1].split(":", 1)
                    hp1 = vals1.strip().split(" / ")
                    hp2 = vals2.strip().split(" / ")
                    h3k4me3 = f"({t1}) {format_encode_value(hp1[0].strip())} / ({t2}) {format_encode_value(hp2[0].strip())}" if len(hp1) >= 1 and len(hp2) >= 1 else "—"
                    h3k27ac = f"({t1}) {format_encode_value(hp1[1].strip())} / ({t2}) {format_encode_value(hp2[1].strip())}" if len(hp1) >= 2 and len(hp2) >= 2 else "—"
                    ctcf = f"({t1}) {format_encode_value(hp1[2].strip())} / ({t2}) {format_encode_value(hp2[2].strip())}" if len(hp1) >= 3 and len(hp2) >= 3 else "—"
                elif len(lines) == 1:
                    if ":" in lines[0]:
                        _, vals = lines[0].split(":", 1)
                    else:
                        vals = lines[0]
                    hp = vals.strip().split(" / ")
                    if len(hp) >= 3:
                        h3k4me3 = format_encode_value(hp[0].strip())
                        h3k27ac = format_encode_value(hp[1].strip())
                        ctcf = format_encode_value(hp[2].strip())
                    elif len(hp) == 2:
                        h3k4me3 = format_encode_value(hp[0].strip())
                        h3k27ac = format_encode_value(hp[1].strip())
            
            # Parse DNase/ATAC values - handle tissue-specific (multi-line) or single-line
            dnase_atac_raw = str(split_gene_value(dnase_atac, g)).strip()
            dnase = "—"
            atac = "—"
            if dnase_atac_raw not in ("nan", "—", ""):
                lines = [l.strip() for l in dnase_atac_raw.split("\n") if l.strip()]
                if len(lines) >= 2:
                    t1 = translate_pt(lines[0].split(":", 1)[0].strip())
                    t2 = translate_pt(lines[1].split(":", 1)[0].strip())
                    _, vals1 = lines[0].split(":", 1)
                    _, vals2 = lines[1].split(":", 1)
                    dp1 = vals1.strip().split(" / ")
                    dp2 = vals2.strip().split(" / ")
                    dnase = f"({t1}) {format_encode_value(dp1[0].strip())} / ({t2}) {format_encode_value(dp2[0].strip())}" if len(dp1) >= 1 and len(dp2) >= 1 else "—"
                    atac = f"({t1}) {format_encode_value(dp1[1].strip())} / ({t2}) {format_encode_value(dp2[1].strip())}" if len(dp1) >= 2 and len(dp2) >= 2 else "—"
                elif len(lines) == 1:
                    if ":" in lines[0]:
                        _, vals = lines[0].split(":", 1)
                    else:
                        vals = lines[0]
                    dp = vals.strip().split(" / ")
                    dnase = format_encode_value(dp[0].strip()) if len(dp) >= 1 else "—"
                    atac = format_encode_value(dp[1].strip()) if len(dp) >= 2 else "—"
            
            records.append({
                "Gene": g,
                "Region": v.get("Region", "—"),
                "Chr": v.get("Chr", "—"),
                "Start": v.get("Start", "—"),
                "Motif": v.get("Motif", "—"),
                "Motif_Ext": motif_ext,
                "Purity": purity,
                "ExpHet": exp_het,
                "RepeatMasker": rm_div,
                "cCRE": ccre,
                "GeneHancer_ID": GENEHANCER_IDS.get(g, "—"),
                "GeneHancer_Score": GENEHANCER_SCORES.get(g, "—"),
                "DNase": dnase,
                "ATAC": atac,
                "H3K4me3": h3k4me3,
                "H3K27ac": h3k27ac,
                "CTCF": ctcf,
                "JARVIS": jarvis_clean,
                "DepletionRank": depletion,
            })
    return pd.DataFrame(records)

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: build_gt
# ─────────────────────────────────────────────────────────────────────────────
def build_gt(df):
    return (
        GT(df, id="tabela_ucsc")
        .tab_header(
            title=html("<strong>Epigenetic and genomic landscape of STRs outliers in genes related to COVID-19</strong>"),
        )
        .tab_spanner(
            label=html("<strong>Genomic Location</strong>"),
            columns=["Gene", "Region", "Chr", "Start", "Motif"],
        )
        .tab_spanner(
            label=html("<strong>TRExplorer</strong>"),
            columns=["Motif_Ext", "Purity", "ExpHet"],
        )
        .tab_spanner(
            label=html("<strong>ENCODE<sup>d</sup></strong>"),
            columns=["DNase", "ATAC", "H3K4me3", "H3K27ac", "CTCF"],
        )
        .cols_label(
            Gene="Gene", Region="Region", Chr="Chr", Start="Start",
            Motif="Motif", Motif_Ext=html("Motif + Ext.<sup>a</sup>"),
            Purity=html("Purity<sup>a</sup>"),
            ExpHet=html("Exp. Het<sup>a</sup>"),
            RepeatMasker=html("RepeatMasker<sup>b</sup>"),
            cCRE=html("cCRE Identity<sup>c</sup>"),
            GeneHancer_ID=html("GeneHancer ID<sup>g</sup>"),
            GeneHancer_Score=html("Score (Value)<sup>g</sup>"),
            DNase="DNase Max",
            ATAC="ATAC Max",
            H3K4me3="H3K4me3",
            H3K27ac="H3K27ac",
            CTCF="CTCF",
            JARVIS=html("JARVIS<sup>e</sup>"),
            DepletionRank=html("DepletionRank<sup>f</sup>"),
        )
        .tab_source_note(
            html("<sup>a</sup> TRExplorer: motif + extension length (repeat copies), repeat purity index "
                 "(fractional sequence perfection; 1.00 = 100% pure), and expected heterozygosity "
                 "derived from the TRExplorer Illumina 174k polymorphic dataset.")
        )
        .tab_source_note(
            html("<sup>b</sup> RepeatMasker: repeat class and divergence from the RepeatMasker annotation.")
        )
        .tab_source_note(
            html("<sup>c</sup> cCRE Identity: candidate cis-Regulatory Element class from ENCODE4 "
                 "(dCRE = distal enhancer).")
        )
        .tab_source_note(
            html("<sup>d</sup> ENCODE signal values are tissue-specific, "
                 "reflecting the tissue in which differential expression "
                 "was identified in RNA-seq analysis.")
        )
        .tab_source_note(
            html("<sup>e</sup> JARVIS: Junk Region Variation Impact Score, reported as Maximum – Minimum "
                 "impact values across tracks overlapping the variant.")
        )
        .tab_source_note(
            html("<sup>f</sup> DepletionRank: UK Biobank / deCODE Genetics depletion rank score "
                 "(Mean ± SD across evaluated bases).")
        )
        .tab_source_note(
            html("<sup>g</sup> GeneHancer: regulatory element–gene interactions (Double Elite) "
                 "from GeneCards. ID = GeneHancer identifier; Score = interaction confidence; "
                 "Value = interaction strength metric.")
        )
        .tab_options(
            table_font_size="13px",
            table_width="100%",
            heading_align="center",
            heading_title_font_size="16px",
            column_labels_font_weight="bold",
            column_labels_border_top_style="solid",
            column_labels_border_bottom_style="solid",
            table_border_top_style="solid",
            table_border_bottom_style="solid",
            table_body_border_bottom_style="solid",
            row_group_padding="6px",
            data_row_padding="5px",
        )
    )

# ─────────────────────────────────────────────────────────────────────────────
# MODULE: save_outputs
# ─────────────────────────────────────────────────────────────────────────────
def save_outputs(df, gt_obj, out_dir, stem="tabela_principal_com_variantes"):
    out_dir = Path(out_dir)
    html_path = out_dir / f"{stem}.html"
    csv_path  = out_dir / f"{stem}.csv"
    xlsx_path = out_dir / f"{stem}.xlsx"

    html_path.write_text(gt_obj.as_raw_html(), encoding="utf-8")
    df.to_csv(csv_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="STR_Loci")
        ws = writer.sheets["STR_Loci"]
        for i, col in enumerate(df.columns, 1):
            ws.column_dimensions[chr(64 + i) if i <= 26 else "A"].width = max(
                df[col].astype(str).str.len().max(), len(col)
            ) + 2

    print(f"HTML → {html_path}")
    print(f"CSV  → {csv_path}")
    print(f"XLSX → {xlsx_path}")

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Build UCSC STR outlier annotation table.")
    parser.add_argument(
        "--outlier-html",
        default=str(BASE / "outlier_detail_table.html"),
        help="Path to the outlier_detail_table.html generated by 6_outlier_detail_table.R",
    )
    args = parser.parse_args()

    SCRIPT_DIR = Path(__file__).resolve().parent
    variant_map = parse_outlier(Path(args.outlier_html))
    df = build_dataframe(BASE / "tabela_principal_estruturada.csv", variant_map)
    gt = build_gt(df)
    save_outputs(df, gt, SCRIPT_DIR)
