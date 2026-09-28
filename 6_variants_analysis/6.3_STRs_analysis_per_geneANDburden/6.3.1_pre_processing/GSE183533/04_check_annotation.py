import pandas as pd
import os

ARQUIVO = "DEG_results/DEG_COVID_vs_CONTROL_FDR_anotado.tsv"

print("=" * 75)
print("CHECAGEM FINAL DOS DEGs - GSE183533")
print("=" * 75)

# 1. LEITURA

print("\nLendo arquivo...")

df = pd.read_csv(ARQUIVO, sep="\t")

print(f"Genes/linhas: {len(df)}")
print(f"Colunas: {len(df.columns)}")

print("\nColunas disponíveis:")
print(df.columns.tolist())


# 2. TOTAL DE DEGs

print("\n" + "=" * 75)
print("1. TOTAL DE DEGs")
print("=" * 75)

print(f"Total: {len(df)}")


# 3. UP / DOWN

print("\n" + "=" * 75)
print("2. DIREÇÃO DA EXPRESSÃO")
print("=" * 75)

if "Direction" in df.columns:

    print(df["Direction"].value_counts(dropna=False))

else:

    print("Coluna Direction não encontrada.")


# 4. BIOTIPOS

print("\n" + "=" * 75)
print("3. BIOTIPOS")
print("=" * 75)

if "biotype" in df.columns:

    biotipo = (
        df["biotype"]
        .fillna("")
        .replace("", "SEM_ANOTACAO")
        .value_counts()
    )

    print(biotipo.to_string())

else:

    print("Coluna biotype não encontrada.")


# 5. PROTEIN CODING

print("\n" + "=" * 75)
print("4. PROTEIN CODING")
print("=" * 75)

if "biotype" in df.columns:

    protein = df[
        df["biotype"].fillna("") == "protein_coding"
    ]

    print(f"Protein-coding: {len(protein)}")


# 6. GENES SEM ANOTAÇÃO

print("\n" + "=" * 75)
print("5. GENES SEM ANOTAÇÃO")
print("=" * 75)

if "gene_symbol" in df.columns:

    sem_symbol = df[
        df["gene_symbol"].isna() |
        (df["gene_symbol"].astype(str).str.strip() == "")
    ]

    print(f"Sem gene symbol: {len(sem_symbol)}")

    if len(sem_symbol) > 0:

        print("\nPrimeiros genes sem símbolo:")

        print(
            sem_symbol[
                ["gene_id", "gene_id_clean", "biotype"]
            ].head(30).to_string(index=False)
        )


# 7. PRINCIPAIS GENES UP

print("\n" + "=" * 75)
print("6. PRINCIPAIS GENES UP")
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


# 8. PRINCIPAIS GENES DOWN

print("\n" + "=" * 75)
print("7. PRINCIPAIS GENES DOWN")
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


# 9. DUPLICIDADE DE GENES

print("\n" + "=" * 75)
print("8. DUPLICIDADE")
print("=" * 75)

if "gene_id_clean" in df.columns:

    duplicados = df["gene_id_clean"].duplicated().sum()

    print(f"Gene IDs duplicados: {duplicados}")

    if duplicados > 0:

        print("\nGenes duplicados:")

        print(
            df[
                df["gene_id_clean"].duplicated(
                    keep=False
                )
            ][
                ["gene_id", "gene_id_clean"]
            ].head(30).to_string(index=False)
        )


# 10. VERIFICAÇÃO DE IDs

print("\n" + "=" * 75)
print("9. VERIFICAÇÃO DOS IDENTIFICADORES")
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
    print(f"Não-ENSG: {len(nao_ensg)}")

    if len(nao_ensg) > 0:

        print("\nIDs não-ENSG:")

        print(
            nao_ensg[
                ["gene_id", "gene_id_clean"]
            ].head(50).to_string(index=False)
        )


# 11. VERIFICAÇÃO DE POSSÍVEIS GENES VIRAIS/MICROBIANOS

print("\n" + "=" * 75)
print("10. VERIFICAÇÃO DE GENES NÃO HUMANOS")
print("=" * 75)

print(
    """
A matriz utilizada contém IDs ENSG (Ensembl Homo sapiens).
Portanto, os resultados desta análise correspondem à
matriz de expressão gênica humana.

A identificação de sequências virais/microbianas deve ser
feita separadamente a partir da análise de reads/metagenômica,
e não misturada com esta matriz de expressão humana.
"""
)


# 12. RESUMO FINAL

print("\n" + "=" * 75)
print("RESUMO FINAL")
print("=" * 75)

print(f"Total de DEGs: {len(df)}")

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
        f"Sem gene symbol: "
        f"{df['gene_symbol'].isna().sum() + (df['gene_symbol'].astype(str).str.strip() == '').sum()}"
    )

if "gene_id_clean" in df.columns:

    print(
        f"IDs ENSG: "
        f"{df['gene_id_clean'].astype(str).str.startswith('ENSG').sum()}"
    )

print("\nChecagem concluída.")
