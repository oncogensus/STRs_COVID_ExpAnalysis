import pandas as pd
import requests
import gzip
import os
import re

# CONFIGURATION

GSE = "GSE183533"

OUTDIR = "metadata"
os.makedirs(OUTDIR, exist_ok=True)

URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE183nnn/GSE183533/matrix/"
    "GSE183533_series_matrix.txt.gz"
)

# DOWNLOAD

print("=" * 75)
print("OFFICIAL METADATA - GSE183533")
print("=" * 75)

arquivo_gz = f"{OUTDIR}/GSE183533_series_matrix.txt.gz"

print("\nDownloading file from GEO...")

r = requests.get(URL, timeout=120)

if r.status_code != 200:
    raise RuntimeError(
        f"Error downloading file: HTTP {r.status_code}"
    )

with open(arquivo_gz, "wb") as f:
    f.write(r.content)

print("Download finished.")

# READ SERIES MATRIX

print("\nReading metadata...")

metadata = {}

with gzip.open(
    arquivo_gz,
    "rt",
    encoding="utf-8",
    errors="replace"
) as f:

    for linha in f:

        linha = linha.rstrip("\n")

        if not linha.startswith("!Sample_"):
            continue

        partes = linha.split("\t")

        campo = partes[0]

        valores = [
            x.strip('"')
            for x in partes[1:]
        ]

        metadata[campo] = valores


# CHECK FIELDS

print("\nFields found:")

for campo in metadata.keys():
    print(" -", campo)

# GSMs

gsm = metadata.get("!Sample_geo_accession", [])

if len(gsm) == 0:
    raise RuntimeError(
        "No GSM found in the Series Matrix."
    )

print(
    f"\nNumber of samples in GEO: {len(gsm)}"
)

# CREATE DATAFRAME

df = pd.DataFrame({
    "GSM": gsm
})

# ADD ALL METADATA

for campo, valores in metadata.items():

    if campo == "!Sample_geo_accession":
        continue

    nome = campo.replace(
        "!Sample_",
        ""
    )

    # ensure same length
    if len(valores) == len(df):

        df[nome] = valores

# IDENTIFY GROUP

print("\nIdentifying groups...")

df["group"] = "UNKNOWN"

for i in range(len(df)):

    texto = " ".join(
        str(x)
        for x in df.iloc[i].tolist()
    ).lower()

    if "covid" in texto:

        df.loc[i, "group"] = "COVID"

    elif (
        "normal" in texto
        or "healthy" in texto
        or "control" in texto
    ):

        df.loc[i, "group"] = "CONTROL"

# IDENTIFY SEX

df["sex"] = "UNKNOWN"

for i in range(len(df)):

    texto = " ".join(
        str(x)
        for x in df.iloc[i].tolist()
    ).lower()

    if re.search(
        r"\bmale\b|\bman\b",
        texto
    ):

        df.loc[i, "sex"] = "Male"

    elif re.search(
        r"\bfemale\b|\bwoman\b",
        texto
    ):

        df.loc[i, "sex"] = "Female"

# SAVE FULL

arquivo_completo = (
    f"{OUTDIR}/GSE183533_metadata_completo.tsv"
)

df.to_csv(
    arquivo_completo,
    sep="\t",
    index=False
)

# SUMMARY

colunas_resumo = [
    "GSM",
    "group",
    "sex"
]

resumo = df[
    [
        c for c in colunas_resumo
        if c in df.columns
    ]
].copy()

arquivo_resumo = (
    f"{OUTDIR}/GSE183533_metadata_resumido.tsv"
)

resumo.to_csv(
    arquivo_resumo,
    sep="\t",
    index=False
)

# RESULTS

print("\n" + "=" * 75)
print("RESULT")
print("=" * 75)

print(
    f"\nSamples found: {len(df)}"
)

print("\nGroups:")

print(
    df["group"]
    .value_counts(dropna=False)
)

print("\nSex:")

print(
    df["sex"]
    .value_counts(dropna=False)
)

print(
    f"\nFull file:"
    f"\n{arquivo_completo}"
)

print(
    f"\nSummary file:"
    f"\n{arquivo_resumo}"
)

print("\nDone.")