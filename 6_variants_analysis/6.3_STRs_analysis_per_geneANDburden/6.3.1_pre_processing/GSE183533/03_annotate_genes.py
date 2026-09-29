import pandas as pd
import requests
import time

ARQUIVO = "DEG_results/DEG_COVID_vs_CONTROL_FDR.tsv"
SAIDA = "DEG_results/DEG_COVID_vs_CONTROL_FDR_anotado.tsv"

print("=" * 70, flush=True)
print("DEG ANNOTATION - GSE183533", flush=True)
print("=" * 70, flush=True)

df = pd.read_csv(ARQUIVO, sep="\t")

print(f"Genes in the file: {len(df)}", flush=True)

df["gene_id_clean"] = (
    df["gene_id"]
    .astype(str)
    .str.replace(r"\..*$", "", regex=True)
)

genes = df["gene_id_clean"].tolist()

SERVER = "https://rest.ensembl.org/lookup/id"

session = requests.Session()

def consultar_gene(gene):

    try:
        resposta = session.get(
            f"{SERVER}/{gene}",
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json"
            },
            timeout=20
        )

        if resposta.status_code == 200:

            info = resposta.json()

            if isinstance(info, dict):
                return {
                    "gene_id_clean": gene,
                    "gene_symbol": info.get("display_name", ""),
                    "biotype": info.get("biotype", ""),
                    "description": info.get("description", "")
                }

    except Exception:
        pass

    return {
        "gene_id_clean": gene,
        "gene_symbol": "",
        "biotype": "",
        "description": ""
    }


resultados = []

TAMANHO_LOTE = 50

for inicio in range(0, len(genes), TAMANHO_LOTE):

    lote = genes[inicio:inicio + TAMANHO_LOTE]

    fim = min(inicio + TAMANHO_LOTE, len(genes))

    print(
        f"Querying genes {inicio + 1}-{fim}/{len(genes)}",
        flush=True
    )

    try:

        resposta = session.post(
            SERVER,
            json={"ids": lote},
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json"
            },
            timeout=30
        )

        if resposta.status_code == 200:

            dados = resposta.json()

            if isinstance(dados, dict):

                for gene in lote:

                    info = dados.get(gene)

                    if isinstance(info, dict):

                        resultados.append({
                            "gene_id_clean": gene,
                            "gene_symbol": info.get("display_name", ""),
                            "biotype": info.get("biotype", ""),
                            "description": info.get("description", "")
                        })

                    else:

                        resultados.append(
                            consultar_gene(gene)
                        )

                continue

    except Exception as e:

        print(
            f"Batch raised an error: {e}",
            flush=True
        )

    print(
        f"Querying {len(lote)} genes individually...",
        flush=True
    )

    for gene in lote:

        resultados.append(
            consultar_gene(gene)
        )

        time.sleep(0.1)


annot = pd.DataFrame(resultados)

df = df.merge(
    annot,
    on="gene_id_clean",
    how="left"
)

df.to_csv(
    SAIDA,
    sep="\t",
    index=False
)

print("\n" + "=" * 70, flush=True)
print("DONE", flush=True)
print("=" * 70, flush=True)

print(f"File: {SAIDA}", flush=True)
print(f"Genes in the result: {len(df)}", flush=True)

print("\nBiotypes:", flush=True)

print(
    df["biotype"]
    .fillna("not_annotated")
    .value_counts()
    .head(20),
    flush=True
)