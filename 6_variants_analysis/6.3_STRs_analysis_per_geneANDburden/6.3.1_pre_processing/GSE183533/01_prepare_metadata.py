import pandas as pd
import requests
import gzip
import os
import re

# CONFIGURAÇÃO

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
print("METADADOS OFICIAIS - GSE183533")
print("=" * 75)

arquivo_gz = f"{OUTDIR}/GSE183533_series_matrix.txt.gz"

print("\nBaixando arquivo do GEO...")

r = requests.get(URL, timeout=120)

if r.status_code != 200:
    raise RuntimeError(
        f"Erro ao baixar arquivo: HTTP {r.status_code}"
    )

with open(arquivo_gz, "wb") as f:
    f.write(r.content)

print("Download concluído.")

# LER SERIES MATRIX

print("\nLendo metadados...")

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


# VERIFICAR CAMPOS

print("\nCampos encontrados:")

for campo in metadata.keys():
    print(" -", campo)

# GSMs

gsm = metadata.get("!Sample_geo_accession", [])

if len(gsm) == 0:
    raise RuntimeError(
        "Nenhum GSM encontrado no Series Matrix."
    )

print(
    f"\nNúmero de amostras no GEO: {len(gsm)}"
)

# CRIAR DATAFRAME

df = pd.DataFrame({
    "GSM": gsm
})

# ADICIONAR TODOS OS METADADOS

for campo, valores in metadata.items():

    if campo == "!Sample_geo_accession":
        continue

    nome = campo.replace(
        "!Sample_",
        ""
    )

    # garantir mesmo tamanho
    if len(valores) == len(df):

        df[nome] = valores

# IDENTIFICAR GROUP

print("\nIdentificando grupos...")

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

# IDENTIFICAR SEXO

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

# SALVAR COMPLETO

arquivo_completo = (
    f"{OUTDIR}/GSE183533_metadata_completo.tsv"
)

df.to_csv(
    arquivo_completo,
    sep="\t",
    index=False
)

# RESUMO

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

# RESULTADOS

print("\n" + "=" * 75)
print("RESULTADO")
print("=" * 75)

print(
    f"\nAmostras encontradas: {len(df)}"
)

print("\nGrupos:")

print(
    df["group"]
    .value_counts(dropna=False)
)

print("\nSexo:")

print(
    df["sex"]
    .value_counts(dropna=False)
)

print(
    f"\nArquivo completo:"
    f"\n{arquivo_completo}"
)

print(
    f"\nArquivo resumido:"
    f"\n{arquivo_resumo}"
)

print("\nConcluído.")
