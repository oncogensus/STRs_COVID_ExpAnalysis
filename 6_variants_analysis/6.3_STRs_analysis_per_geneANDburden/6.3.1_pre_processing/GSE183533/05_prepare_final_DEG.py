import pandas as pd
import os

arquivo = "DEG_results/DEG_COVID_vs_CONTROL_FDR_anotado.tsv"
saida = "DEG_results"

df = pd.read_csv(arquivo, sep="\t")

# 1. FULL FILE

df.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_FINAL.tsv",
    sep="\t",
    index=False
)

# 2. PROTEIN-CODING

protein = df[
    df["biotype"].fillna("") == "protein_coding"
].copy()

protein.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_protein_coding.tsv",
    sep="\t",
    index=False
)

# 3. UP

up = df[
    df["Direction"] == "Up_COVID"
].copy()

up.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_UP_FINAL.tsv",
    sep="\t",
    index=False
)

# 4. DOWN

down = df[
    df["Direction"] == "Down_COVID"
].copy()

down.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_DOWN_FINAL.tsv",
    sep="\t",
    index=False
)

# 5. PROTEIN-CODING UP

protein_up = protein[
    protein["Direction"] == "Up_COVID"
].copy()

protein_up.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_protein_coding_UP.tsv",
    sep="\t",
    index=False
)

# 6. PROTEIN-CODING DOWN

protein_down = protein[
    protein["Direction"] == "Down_COVID"
].copy()

protein_down.to_csv(
    f"{saida}/GSE183533_DEG_COVID_vs_CONTROL_protein_coding_DOWN.tsv",
    sep="\t",
    index=False
)

# SUMMARY

print("=" * 70)
print("FINAL FILES - GSE183533")
print("=" * 70)

print(f"Total DEGs:            {len(df)}")
print(f"UP COVID:              {len(up)}")
print(f"DOWN COVID:            {len(down)}")
print(f"Protein-coding:        {len(protein)}")
print(f"Protein-coding UP:     {len(protein_up)}")
print(f"Protein-coding DOWN:   {len(protein_down)}")

print("\nGenerated files:")

for arquivo in sorted(os.listdir(saida)):
    if "GSE183533" in arquivo:
        print(" -", arquivo)

print("\nDone.")