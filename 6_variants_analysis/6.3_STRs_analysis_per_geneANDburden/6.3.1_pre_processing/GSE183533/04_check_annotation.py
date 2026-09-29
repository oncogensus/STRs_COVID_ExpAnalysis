import pandas as pd
import os

ARQUIVO = "DEG_results/DEG_COVID_vs_CONTROL_FDR_anotado.tsv"

print("=" * 75)
print("FINAL DEG CHECK - GSE183533")
print("=" * 75)

# 1. READ

print("\nReading file...")

df = pd.read_csv(ARQUIVO, sep="\t")

print(f"Genes/rows: {len(df)}")
print(f"Columns: {len(df.columns)}")

print("\nAvailable columns:")
print(df.columns.tolist())


# 2. TOTAL DEGs

print("\n" + "=" * 75)
print("1. TOTAL DEGs")
print("=" * 75)

print(f"Total: {len(df)}")


# 3. UP / DOWN

print("\n" + "=" * 75)
print("2. EXPRESSION DIRECTION")
print("=" * 75)

if "Direction" in df.columns:

    print(df["Direction"].value_counts(dropna=False))

else:

    print("Column Direction not found.")


# 4. BIOTYPES

print("\n" + "=" * 75)
print("3. BIOTYPES")
print("=" * 75)

if "biotype" in df.columns:

    biotipo = (
        df["biotype"]
        .fillna("")
        .replace("", "NO_ANNOTATION")
        .value_counts()
    )

    print(biotipo.to_string())

else:

    print("Column biotype not found.")


# 5. PROTEIN CODING

print("\n" + "=" * 75)
print("4. PROTEIN CODING")
print("=" * 75)

if "biotype" in df.columns:

    protein = df[
        df["biotype"].fillna("") == "protein_coding"
    ]

    print(f"Protein-coding: {len(protein)}")


# 6. GENES WITHOUT ANNOTATION

print("\n" + "=" * 75)
print("5. GENES WITHOUT ANNOTATION")
print("=" * 75)

if "gene_symbol" in df.columns:

    sem_symbol = df[
        df["gene_symbol"].isna() |
        (df["gene_symbol"].astype(str).str.strip() == "")
    ]

    print(f"Without gene symbol: {len(sem_symbol)}")

    if len(sem_symbol) > 0:

        print("\nFirst genes without symbol:")

        print(
            sem_symbol[
                ["gene_id", "gene_id_clean", "biotype"]
            ].head(30).to_string(index=False)
        )


# 7. TOP UP GENES

print("\n" + "=" * 75)
print("6. TOP UP GENES")
print("=" * 75)

if "Direction" in df.columns:

    up = df[
        df["Direction"].astype(str).str.contains(
            "UP",
            case=False,
            na=False
        )
    ].copy()

    print(f"Total UP: {len(up)}")

    if "FDR" in up.columns:

        up = up.sort_values("FDR")

    colunas = [
        c for c in
        ["gene_symbol", "gene_id", "logFC", "FDR", "biotype"]
        if c in up.columns
    ]

    print(
        up[colunas].head(30).to_string(index=False)
    )


# 8. TOP DOWN GENES

print("\n" + "=" * 75)
print("7. TOP DOWN GENES")
print("=" * 75)

if "Direction" in df.columns:

    down = df[
        df["Direction"].astype(str).str.contains(
            "DOWN",
            case=False,
            na=False
        )
    ].copy()

    print(f"Total DOWN: {len(down)}")

    if "FDR" in down.columns:

        down = down.sort_values("FDR")

    colunas = [
        c for c in
        ["gene_symbol", "gene_id", "logFC", "FDR", "biotype"]
        if c in down.columns
    ]

    print(
        down[colunas].head(30).to_string(index=False)
    )


# 9. GENE DUPLICATES

print("\n" + "=" * 75)
print("8. DUPLICATES")
print("=" * 75)

if "gene_id_clean" in df.columns:

    duplicados = df["gene_id_clean"].duplicated().sum()

    print(f"Duplicated gene IDs: {duplicados}")

    if duplicados > 0:

        print("\nDuplicated genes:")

        print(
            df[
                df["gene_id_clean"].duplicated(
                    keep=False
                )
            ][
                ["gene_id", "gene_id_clean"]
            ].head(30).to_string(index=False)
        )


# 10. ID VERIFICATION

print("\n" + "=" * 75)
print("9. IDENTIFIER VERIFICATION")
print("=" * 75)

if "gene_id_clean" in df.columns:

    ensg = df[
        df["gene_id_clean"]
        .astype(str)
        .str.startswith("ENSG")
    ]

    nao_ensg = df[
        ~df["gene_id_clean"]
        .astype(str)
        .str.startswith("ENSG")
    ]

    print(f"ENSG: {len(ensg)}")
    print(f"Non-ENSG: {len(nao_ensg)}")

    if len(nao_ensg) > 0:

        print("\nNon-ENSG IDs:")

        print(
            nao_ensg[
                ["gene_id", "gene_id_clean"]
            ].head(50).to_string(index=False)
        )


# 11. CHECK FOR POSSIBLE VIRAL/MICROBIAL GENES

print("\n" + "=" * 75)
print("10. NON-HUMAN GENE VERIFICATION")
print("=" * 75)

print(
    """
The matrix used contains ENSG IDs (Ensembl Homo sapiens).
Therefore, the results of this analysis correspond to the
human gene expression matrix.

The identification of viral/microbial sequences must be
done separately from the analysis of reads/metagenomics,
and should not be mixed with this human expression matrix.
"""
)


# 12. FINAL SUMMARY

print("\n" + "=" * 75)
print("FINAL SUMMARY")
print("=" * 75)

print(f"Total DEGs: {len(df)}")

if "Direction" in df.columns:

    print(
        f"UP: {(df['Direction'].astype(str).str.contains('UP', case=False)).sum()}"
    )

    print(
        f"DOWN: {(df['Direction'].astype(str).str.contains('DOWN', case=False)).sum()}"
    )

if "biotype" in df.columns:

    print(
        f"Protein-coding: "
        f"{(df['biotype'].fillna('') == 'protein_coding').sum()}"
    )

if "gene_symbol" in df.columns:

    print(
        f"Without gene symbol: "
        f"{df['gene_symbol'].isna().sum() + (df['gene_symbol'].astype(str).str.strip() == '').sum()}"
    )

if "gene_id_clean" in df.columns:

    print(
        f"ENSG IDs: "
        f"{df['gene_id_clean'].astype(str).str.startswith('ENSG').sum()}"
    )

print("\nCheck finished.")