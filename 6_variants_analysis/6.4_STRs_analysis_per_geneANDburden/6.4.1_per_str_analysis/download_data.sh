#!/usr/bin/env bash
# download_data.sh
# Baixa os summary stats completos do COVID-19 Host Genetics Initiative r7
# (freeze 7, 2022-04-03; arquivos "leave_23andme" = completos, sem 23andMe)
# e a anotacao de genes GENCODE v39 (hg38) para mapear variantes->genes.
#
# Rodar no cluster (carloschagas), no micromamba env 'igv' (tem python base).
#   bash download_data.sh
# Os arquivos vao para ./data (mesmo dir deste script).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/data"
mkdir -p "$DATA"

# ---- COVID-19 HG r7 full summary stats (A2 infeccao, B2 hospitalizacao, C2 critico) ----
BASE="https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats"
FILES=(
  "COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz"
  "COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz"
  "COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz"
)

echo "=== Baixando COVID-19 HG r7 summary stats (leave_23andme) ==="
for f in "${FILES[@]}"; do
  out="$DATA/$f"
  if [ -s "$out" ]; then
    echo "  [ja existe] $f"
  else
    echo "  [baixando] $f"
    curl -L -C - -o "$out" "$BASE/$f" || echo "  AVISO: falha ao baixar $f"
  fi
done

# ---- GENCODE v39 primary_assembly GTF (hg38, com prefixo chr) ----
GTF="gencode.v39.primary_assembly.annotation.gtf.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/$GTF"
if [ -s "$DATA/$GTF" ]; then
  echo "[ja existe] $GTF"
else
  echo "[baixando] $GTF"
  curl -L -C - -o "$DATA/$GTF" "$GTF_URL" || echo "  AVISO: falha ao baixar $GTF"
fi

echo
echo "=== Arquivos em $DATA ==="
ls -lh "$DATA"
echo
echo "Se algum download falhou (sem internet no cluster), baixe no seu PC e rsync/scp para $DATA."
echo "Build: rode depois  bash build_gene_bed.py  e  bash extract_covid_genes.py"
