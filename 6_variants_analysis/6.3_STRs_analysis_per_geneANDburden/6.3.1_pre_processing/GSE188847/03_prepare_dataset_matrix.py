#!/usr/bin/env python3

import os
import gzip
import re
import pandas as pd

SALMON_DIR = "salmon"
OUT_MATRIX = "GSE188847_matrix.tsv.gz"
OUT_METADATA = "GSE188847_metadata.tsv"

# Selecionar somente COVID, CONTROL e ICUVENT
arquivos = []

for f in sorted(os.listdir(SALMON_DIR)):

    if not f.endswith(".gz"):
        continue

    # Somente os 54 pacientes desejados
    m = re.match(r"(GSM[0-9]+)_(COVID|CONTROL|ICUVENT)([0-9]+)\.salmon\.", f)

    if m:
        gsm = m.group(1)
        grupo = m.group(2)
        numero = m.group(3)

        arquivos.append({
            "file": f,
            "GSM": gsm,
            "sample": grupo + numero,
            "group": grupo
        })

print("=" * 70)
print("GSE188847")
print("=" * 70)

print("Arquivos selecionados:", len(arquivos))

for grupo in ["COVID", "CONTROL", "ICUVENT"]:
    n = sum(x["group"] == grupo for x in arquivos)
    print("{}: {}".format(grupo, n))

print("=" * 70)

if len(arquivos) != 54:
    raise Exception(
        "ERRO: esperadas 54 amostras, mas foram encontradas {}".format(
            len(arquivos)
        )
    )

# ---------------------------------------------------------
# Ler arquivos Salmon
# ---------------------------------------------------------

matriz = None

for i, info in enumerate(arquivos):

    caminho = os.path.join(SALMON_DIR, info["file"])

    print(
        "[{}/{}] {}".format(
            i + 1,
            len(arquivos),
            info["file"]
        )
    )

    with gzip.open(caminho, "rt") as f:

        # Detectar automaticamente o formato
        df = pd.read_csv(f, sep="\t")

    # Salmon pode ter:
    # Name / Length / EffectiveLength / TPM / NumReads
    if "Name" not in df.columns:
        raise Exception(
            "Coluna Name não encontrada em {}".format(info["file"])
        )

    if "NumReads" not in df.columns:
        raise Exception(
            "Coluna NumReads não encontrada em {}".format(info["file"])
        )

    # Manter apenas gene/transcrito e contagem
    df = df[["Name", "NumReads"]].copy()

    df.columns = ["gene", info["sample"]]

    # Remover duplicatas
    df = df.groupby("gene", as_index=False)[info["sample"]].sum()

    if matriz is None:
        matriz = df
    else:
        matriz = matriz.merge(
            df,
            on="gene",
            how="outer"
        )

# ---------------------------------------------------------
# Limpeza
# ---------------------------------------------------------

matriz = matriz.fillna(0)

# Ordenar genes
matriz = matriz.sort_values("gene")

# ---------------------------------------------------------
# Salvar matriz
# ---------------------------------------------------------

print("=" * 70)
print("Salvando matriz...")

matriz.to_csv(
    OUT_MATRIX,
    sep="\t",
    index=False,
    compression="gzip"
)

# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------

metadata = pd.DataFrame(arquivos)

metadata = metadata[
    ["GSM", "sample", "group", "file"]
]

metadata.to_csv(
    OUT_METADATA,
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# Resumo
# ---------------------------------------------------------

print("=" * 70)
print("CONCLUÍDO")
print("=" * 70)

print("Matriz:", OUT_MATRIX)
print("Metadata:", OUT_METADATA)

print()
print("Dimensão da matriz:")
print("Genes:", matriz.shape[0])
print("Amostras:", matriz.shape[1] - 1)

print()
print("Grupos:")
print(metadata["group"].value_counts())

print("=" * 70)
