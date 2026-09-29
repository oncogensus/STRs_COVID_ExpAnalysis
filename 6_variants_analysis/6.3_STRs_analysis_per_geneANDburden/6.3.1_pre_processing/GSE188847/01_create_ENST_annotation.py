import gzip
import re

gtf = "Homo_sapiens.GRCh38.104.gtf.gz"
out = "ENST_to_gene.tsv"

seen = set()

with gzip.open(gtf, "rt") as f, open(out, "w") as o:

    o.write("transcript_id\tgene_id\tgene_name\n")

    for line in f:

        if line.startswith("#"):
            continue

        if '\ttranscript\t' not in line:
            continue

        fields = line.rstrip("\n").split("\t")
        attr = fields[8]

        transcript = re.search(r'transcript_id "([^"]+)"', attr)
        gene_id = re.search(r'gene_id "([^"]+)"', attr)
        gene_name = re.search(r'gene_name "([^"]+)"', attr)

        if transcript and gene_id:

            tid = transcript.group(1).split(".")[0]
            gid = gene_id.group(1).split(".")[0]

            if gene_name:
                gname = gene_name.group(1)
            else:
                gname = gid

            if tid not in seen:
                o.write(tid + "\t" + gid + "\t" + gname + "\n")
                seen.add(tid)

print("Done.")
print("Annotated transcripts:", len(seen))