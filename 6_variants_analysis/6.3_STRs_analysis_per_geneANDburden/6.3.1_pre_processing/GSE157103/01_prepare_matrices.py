#!/usr/bin/env python3

import gzip
import csv
import re
import sys

# FILES

MATRIX = "GSE157103_genes.ec.tsv.gz"
METADATA = "GSE157103_clinical_metadata.tsv"

OUT_COVID_MATRIX = "COVID_ICU_NonICU_matrix.tsv.gz"
OUT_NONCOVID_MATRIX = "NONCOVID_ICU_NonICU_matrix.tsv.gz"

OUT_COVID_META = "COVID_ICU_NonICU_metadata.tsv"
OUT_NONCOVID_META = "NONCOVID_ICU_NonICU_metadata.tsv"


# 1. READ METADATA

print("\n==============================================")
print("GSE157103 - MATRIX PREPARATION")
print("==============================================")

metadata = []

with open(METADATA, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata.append(row)

print("\nTotal samples in metadata:", len(metadata))


# 2. SPLIT COVID / NONCOVID

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


# 3. FUNCTION TO EXTRACT NUMBER FROM SAMPLE TITLE

def extract_number(title):

    # COVID_77
    # NONCOVID_22

    m = re.search(r"_(\d+)_", title)

    if not m:
        raise ValueError(
            "Could not extract number from: {}".format(title)
        )

    return int(m.group(1))


# 4. CREATE MAPPING

covid_map = {}

for row in covid:

    n = extract_number(row["sample_title"])

    key = "C{}".format(n)

    if key in covid_map:
        print("ERROR: duplicate:", key)
        sys.exit(1)

    covid_map[key] = row


noncovid_map = {}

for row in noncovid:

    n = extract_number(row["sample_title"])

    key = "NC{}".format(n)

    if key in noncovid_map:
        print("ERROR: duplicate:", key)
        sys.exit(1)

    noncovid_map[key] = row


# 5. VERIFY NUMBERING

print("\n==============================================")
print("VERIFYING MAPPING")
print("==============================================")


print("\nCOVID examples:")

for key in sorted(covid_map.keys(), key=lambda x: int(x[1:]))[:15]:

    row = covid_map[key]

    print(
        "{} -> {} -> {}".format(
            key,
            row["GSM"],
            row["sample_title"]
        )
    )


print("\nNONCOVID examples:")

for key in sorted(noncovid_map.keys(), key=lambda x: int(x[2:]))[:15]:

    row = noncovid_map[key]

    print(
        "{} -> {} -> {}".format(
            key,
            row["GSM"],
            row["sample_title"]
        )
    )


# 6. VERIFY MATRIX

print("\n==============================================")
print("VERIFYING EXPRESSION MATRIX")
print("==============================================")

with gzip.open(MATRIX, "rt") as f:

    reader = csv.reader(f, delimiter="\t")

    header = next(reader)

    original_samples = header[1:]

    print("Number of samples in matrix:", len(original_samples))

    if len(original_samples) != 126:

        print("ERROR: expected 126 samples.")
        sys.exit(1)


    # Verify that all C/NC exist

    missing = []

    for key in covid_map:

        if key not in original_samples:
            missing.append(key)

    for key in noncovid_map:

        if key not in original_samples:
            missing.append(key)


    if missing:

        print("\nERROR: samples not found in matrix:")

        for x in sorted(missing):
            print(x)

        sys.exit(1)


    print("All identifiers were found.")


    # ========================================================
    # 7. DEFINE SAMPLE ORDER
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
    # 8. CONVERT C/NC TO MATRIX INDICES
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
    # 9. GSMs IN THE SAME ORDER AS THE NEW MATRIX
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
    # 10. SAVE MATRICES
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
    # 11. PROCESS GENES
    # ========================================================

    genes = 0

    for row in reader:

        if len(row) != 127:

            print(
                "\nERROR: line with an incorrect number of columns."
            )

            print(
                "Gene:",
                row[0]
            )

            print(
                "Number of fields:",
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
                "Genes processed:",
                genes
            )


    covid_out.close()
    noncovid_out.close()


# 12. SAVE METADATA

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


# 13. FINAL SUMMARY

print("\n==============================================")
print("DONE")
print("==============================================")

print("\nGenes:", genes)

print("\nCOVID:")
print("  Samples:", len(covid_gsms))
print("  Matrix:", OUT_COVID_MATRIX)
print("  Metadata:", OUT_COVID_META)

print("\nNONCOVID:")
print("  Samples:", len(noncovid_gsms))
print("  Matrix:", OUT_NONCOVID_MATRIX)
print("  Metadata:", OUT_NONCOVID_META)

print("\n==============================================")