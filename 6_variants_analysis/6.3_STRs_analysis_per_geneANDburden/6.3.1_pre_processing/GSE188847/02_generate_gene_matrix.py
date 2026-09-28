import gzip
import os
import re
from collections import defaultdict

SALMON_DIR = "salmon"
ANNOT = "ENST_to_gene.tsv"

OUT = "GSE188847_gene_counts.tsv"
META = "GSE188847_gene_metadata.tsv"

# ==========================================================
# 1. LER ANOTAÇÃO
# ==========================================================

annotation = {}
gene_names = {}

with open(ANNOT, "r") as f:

    next(f)

    for line in f:

        line = line.rstrip("\n")

        if not line:
            continue

        parts = line.split("\t")

        if len(parts) < 3:
            continue

        transcript = parts[0]
        gene_id = parts[1]
        gene_name = parts[2]

        annotation[transcript] = gene_id

        if gene_id not in gene_names:
            gene_names[gene_id] = gene_name

print("=" * 70)
print("ANOTAÇÃO")
print("=" * 70)

print("Transcritos anotados:", len(annotation))
print("Genes:", len(gene_names))


# ==========================================================
# 2. SELECIONAR AMOSTRAS CLÍNICAS
# ==========================================================

files = []

for filename in sorted(os.listdir(SALMON_DIR)):

    if not filename.endswith(".gz"):
        continue

    if re.search(r"_COVID[0-9]+\.salmon\.", filename):

        group = "COVID"

    elif re.search(r"_CONTROL[0-9]+\.salmon\.", filename):

        group = "CONTROL"

    elif re.search(r"_ICUVENT[0-9]+\.salmon\.", filename):

        group = "ICUVENT"

    else:

        continue

    files.append((filename, group))


# Ordenar pelo nome da amostra
def sample_number(item):

    filename = item[0]

    match = re.search(
        r"_(COVID|CONTROL|ICUVENT)([0-9]+)",
        filename
    )

    if match:

        group = match.group(1)
        number = int(match.group(2))

        return (group, number)

    return (filename, 0)


files = sorted(files, key=sample_number)


print("\n" + "=" * 70)
print("AMOSTRAS")
print("=" * 70)

print("Total:", len(files))


# ==========================================================
# 3. MATRIZ
# ==========================================================

matrix = defaultdict(lambda: defaultdict(float))

sample_names = []
sample_groups = {}

for i, (filename, group) in enumerate(files, 1):

    match = re.search(
        r"_(COVID|CONTROL|ICUVENT)([0-9]+)",
        filename
    )

    if not match:
        continue

    sample = match.group(1) + match.group(2)

    sample_names.append(sample)
    sample_groups[sample] = group

    print("[%d/%d] %s" %
          (i, len(files), sample))

    path = os.path.join(SALMON_DIR, filename)

    with gzip.open(path, "rt") as f:

        header = f.readline().rstrip("\n").split("\t")

        idx_name = header.index("Name")
        idx_reads = header.index("NumReads")

        for line in f:

            parts = line.rstrip("\n").split("\t")

            if len(parts) <= max(idx_name, idx_reads):
                continue

            transcript = parts[idx_name]

            # Remover versão do ENST
            transcript = transcript.split(".")[0]

            if transcript not in annotation:
                continue

            gene_id = annotation[transcript]

            try:

                reads = float(parts[idx_reads])

            except:

                continue

            matrix[gene_id][sample] += reads


# ==========================================================
# 4. SALVAR MATRIZ
# ==========================================================

print("\n" + "=" * 70)
print("SALVANDO MATRIZ")
print("=" * 70)

with open(OUT, "w") as o:

    # Cabeçalho
    header = ["gene_id", "gene_name"] + sample_names

    o.write("\t".join(header) + "\n")

    # Genes
    for gene_id in sorted(matrix.keys()):

        gene_name = gene_names.get(
            gene_id,
            gene_id
        )

        values = [
            gene_id,
            gene_name
        ]

        for sample in sample_names:

            value = matrix[gene_id].get(
                sample,
                0
            )

            values.append(
                "%.6f" % value
            )

        o.write(
            "\t".join(values) + "\n"
        )


# ==========================================================
# 5. METADATA
# ==========================================================

with open(META, "w") as o:

    o.write("sample\tgroup\n")

    for sample in sample_names:

        o.write(
            sample + "\t" +
            sample_groups[sample] + "\n"
        )


# ==========================================================
# 6. RESUMO
# ==========================================================

print("\n" + "=" * 70)
print("CONCLUÍDO")
print("=" * 70)

print("Matriz:", OUT)
print("Metadata:", META)
print("Genes:", len(matrix))
print("Amostras:", len(sample_names))

print("\nGrupos:")

counts = {}

for sample in sample_names:

    group = sample_groups[sample]

    if group not in counts:
        counts[group] = 0

    counts[group] += 1

for group in ["COVID", "CONTROL", "ICUVENT"]:

    print(
        "%s: %d" %
        (group, counts.get(group, 0))
    )
