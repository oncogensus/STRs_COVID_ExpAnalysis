#!/usr/bin/env python3

import os
import gzip
import re
import pandas as pd

SALMON_DIR = "salmon"
OUT_MATRIX = "GSE188847_matrix.tsv.gz"
OUT_METADATA = "GSE188847_metadata.tsv"

# Select only COVID, CONTROL and ICUVENT
arquivos = []

for f in sorted(os.listdir(SALMON_DIR)):

    if not f.endswith(".gz"):
        continue

    # Only the 54 desired patients
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

print("Selected files:", len(arquivos))

for grupo in ["COVID", "CONTROL", "ICUVENT"]:
    n = sum(x["group"] == grupo for x in arquivos)
    print("{}: {}".format(grupo, n))

print("=" * 70)

if len(arquivos) != 54:
    raise Exception(
        "ERROR: expected 54 samples, but {} were found".format(
            len(arquivos)
        )
    )

# ---------------------------------------------------------
# Read Salmon files
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

        # Automatically detect the format
        df = pd.read_csv(f, sep="\t")

    # Salmon can have:
    # Name / Length / EffectiveLength / TPM / NumReads
    if "Name" not in df.columns:
        raise Exception(
            "Column Name not found in {}".format(info["file"])
        )

    if "NumReads" not in df.columns:
        raise Exception(
            "Column NumReads not found in {}".format(info["file"])
        )

    # Keep only gene/transcript and count
    df = df[["Name", "NumReads"]].copy()

    df.columns = ["gene", info["sample"]]

    # Remove duplicates
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
# Cleaning
# ---------------------------------------------------------

matriz = matriz.fillna(0)

# Sort genes
matriz = matriz.sort_values("gene")

# ---------------------------------------------------------
# Save matrix
# ---------------------------------------------------------

print("=" * 70)
print("Saving matrix...")

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
# Summary
# ---------------------------------------------------------

print("=" * 70)
print("DONE")
print("=" * 70)

print("Matrix:", OUT_MATRIX)
print("Metadata:", OUT_METADATA)

print()
print("Matrix dimensions:")
print("Genes:", matriz.shape[0])
print("Samples:", matriz.shape[1] - 1)

print()
print("Groups:")
print(metadata["group"].value_counts())

print("=" * 70)