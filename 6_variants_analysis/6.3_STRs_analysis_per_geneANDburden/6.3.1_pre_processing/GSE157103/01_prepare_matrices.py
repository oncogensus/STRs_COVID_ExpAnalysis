#!/usr/bin/env python3

import gzip
import csv
import re
import sys

# ARQUIVOS

MATRIX = "GSE157103_genes.ec.tsv.gz"
METADATA = "GSE157103_clinical_metadata.tsv"

OUT_COVID_MATRIX = "COVID_ICU_NonICU_matrix.tsv.gz"
OUT_NONCOVID_MATRIX = "NONCOVID_ICU_NonICU_matrix.tsv.gz"

OUT_COVID_META = "COVID_ICU_NonICU_metadata.tsv"
OUT_NONCOVID_META = "NONCOVID_ICU_NonICU_metadata.tsv"


# 1. LER METADATA

print("\n==============================================")
print("GSE157103 - PREPARACAO DAS MATRIZES")
print("==============================================")

metadata = []

with open(METADATA, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata.append(row)

print("\nTotal de amostras na metadata:", len(metadata))


# 2. SEPARAR COVID / NONCOVID

covid = [
    x for x in metadata
    if x["condition"] == "COVID"
]

noncovid = [
    x for x in metadata
    if x["condition"] == "NONCOVID"
]


print("\nCOVID")
print("  Total:", len(covid))
print("  ICU:", sum(x["ICU_status"] == "ICU" for x in covid))
print("  NonICU:", sum(x["ICU_status"] == "NonICU" for x in covid))

print("\nNONCOVID")
print("  Total:", len(noncovid))
print("  ICU:", sum(x["ICU_status"] == "ICU" for x in noncovid))
print("  NonICU:", sum(x["ICU_status"] == "NonICU" for x in noncovid))


# 3. FUNCAO PARA EXTRAIR NUMERO DO SAMPLE TITLE

def extract_number(title):

    # COVID_77
    # NONCOVID_22

    m = re.search(r"_(\d+)_", title)

    if not m:
        raise ValueError(
            "Nao foi possivel extrair numero de: {}".format(title)
        )

    return int(m.group(1))


# 4. CRIAR MAPEAMENTO

covid_map = {}

for row in covid:

    n = extract_number(row["sample_title"])

    key = "C{}".format(n)

    if key in covid_map:
        print("ERRO: duplicacao:", key)
        sys.exit(1)

    covid_map[key] = row


noncovid_map = {}

for row in noncovid:

    n = extract_number(row["sample_title"])

    key = "NC{}".format(n)

    if key in noncovid_map:
        print("ERRO: duplicacao:", key)
        sys.exit(1)

    noncovid_map[key] = row


# 5. VERIFICAR NUMERACAO

print("\n==============================================")
print("VERIFICANDO MAPEAMENTO")
print("==============================================")


print("\nExemplos COVID:")

for key in sorted(covid_map.keys(), key=lambda x: int(x[1:]))[:15]:

    row = covid_map[key]

    print(
        "{} -> {} -> {}".format(
            key,
            row["GSM"],
            row["sample_title"]
        )
    )


print("\nExemplos NONCOVID:")

for key in sorted(noncovid_map.keys(), key=lambda x: int(x[2:]))[:15]:

    row = noncovid_map[key]

    print(
        "{} -> {} -> {}".format(
            key,
            row["GSM"],
            row["sample_title"]
        )
    )


# 6. VERIFICAR MATRIZ

print("\n==============================================")
print("VERIFICANDO MATRIZ DE EXPRESSAO")
print("==============================================")

with gzip.open(MATRIX, "rt") as f:

    reader = csv.reader(f, delimiter="\t")

    header = next(reader)

    original_samples = header[1:]

    print("Numero de amostras na matriz:", len(original_samples))

    if len(original_samples) != 126:

        print("ERRO: esperadas 126 amostras.")
        sys.exit(1)


    # Verificar se todos os C/NC existem

    missing = []

    for key in covid_map:

        if key not in original_samples:
            missing.append(key)

    for key in noncovid_map:

        if key not in original_samples:
            missing.append(key)


    if missing:

        print("\nERRO: amostras nao encontradas na matriz:")

        for x in sorted(missing):
            print(x)

        sys.exit(1)


    print("Todos os identificadores foram encontrados.")


    # ========================================================
    # 7. DEFINIR ORDEM DAS AMOSTRAS
    # ========================================================

    covid_keys = sorted(
        covid_map.keys(),
        key=lambda x: int(x[1:])
    )

    noncovid_keys = sorted(
        noncovid_map.keys(),
        key=lambda x: int(x[2:])
    )


    # ========================================================
    # 8. CONVERTER C/NC PARA INDICES DA MATRIZ
    # ========================================================

    index = {
        sample: i
        for i, sample in enumerate(original_samples, start=1)
    }


    covid_indices = [
        index[x]
        for x in covid_keys
    ]

    noncovid_indices = [
        index[x]
        for x in noncovid_keys
    ]


    # ========================================================
    # 9. GSMs NA MESMA ORDEM DA MATRIZ NOVA
    # ========================================================

    covid_gsms = [
        covid_map[x]["GSM"]
        for x in covid_keys
    ]

    noncovid_gsms = [
        noncovid_map[x]["GSM"]
        for x in noncovid_keys
    ]


    # ========================================================
    # 10. SALVAR MATRIZES
    # ========================================================

    covid_out = gzip.open(
        OUT_COVID_MATRIX,
        "wt"
    )

    noncovid_out = gzip.open(
        OUT_NONCOVID_MATRIX,
        "wt"
    )


    covid_out.write(
        "#symbol\t" +
        "\t".join(covid_gsms) +
        "\n"
    )

    noncovid_out.write(
        "#symbol\t" +
        "\t".join(noncovid_gsms) +
        "\n"
    )


    # ========================================================
    # 11. PROCESSAR GENES
    # ========================================================

    genes = 0

    for row in reader:

        if len(row) != 127:

            print(
                "\nERRO: linha com numero incorreto de colunas."
            )

            print(
                "Gene:",
                row[0]
            )

            print(
                "Numero de campos:",
                len(row)
            )

            sys.exit(1)


        gene = row[0]

        covid_values = [
            row[i]
            for i in covid_indices
        ]

        noncovid_values = [
            row[i]
            for i in noncovid_indices
        ]


        covid_out.write(
            gene + "\t" +
            "\t".join(covid_values) +
            "\n"
        )

        noncovid_out.write(
            gene + "\t" +
            "\t".join(noncovid_values) +
            "\n"
        )


        genes += 1

        if genes % 5000 == 0:

            print(
                "Genes processados:",
                genes
            )


    covid_out.close()
    noncovid_out.close()


# 12. SALVAR METADATA

fields = [
    "GSM",
    "sample_title",
    "condition",
    "age",
    "sex",
    "ICU_status"
]


with open(OUT_COVID_META, "w") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for key in covid_keys:

        writer.writerow(
            covid_map[key]
        )


with open(OUT_NONCOVID_META, "w") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for key in noncovid_keys:

        writer.writerow(
            noncovid_map[key]
        )


# 13. RESUMO FINAL

print("\n==============================================")
print("CONCLUIDO")
print("==============================================")

print("\nGenes:", genes)

print("\nCOVID:")
print("  Amostras:", len(covid_gsms))
print("  Matriz:", OUT_COVID_MATRIX)
print("  Metadata:", OUT_COVID_META)

print("\nNONCOVID:")
print("  Amostras:", len(noncovid_gsms))
print("  Matriz:", OUT_NONCOVID_MATRIX)
print("  Metadata:", OUT_NONCOVID_META)

print("\n==============================================")
