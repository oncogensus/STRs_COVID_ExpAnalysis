#!/usr/bin/env bash
# run_igv_option1.sh
# Snapshots IGV (headless) das 8 regioes STR em DUAS tracks:
#   str_samples_with_variant.bed    (amostra COM variante)
#   str_samples_without_variant.bed (controle SEM variante)
# Requer: env micromamba 'igv' (bioconda igv + conda-forge xorg-xserver-xvfb)
#         e acesso a hg38 pela internet (IGV baixa o genoma).
# Uso: bash run_igv_option1.sh   (roda no dir 7.4_scRNA-seq_analysis/)
set -euo pipefail
cd "$(dirname "$0")"

BED_WITH=str_samples_with_variant.bed
BED_WITHOUT=str_samples_without_variant.bed
EXT_DIR=external_tracks
OUTDIR=snapshots_option1

for f in "$BED_WITH" "$BED_WITHOUT"; do
  [ -f "$f" ] || { echo "ERRO: $f ausente. Rode str_samples_to_bed.R primeiro."; exit 1; }
done
mkdir -p "$OUTDIR"

# --- Xvfb (display virtual); reusa DISPLAY se ja houver ---
if [ -z "${DISPLAY:-}" ]; then
  XDISP=:99
  Xvfb "$XDISP" -screen 0 1024x768x24 >/dev/null 2>&1 &
  XVFB_PID=$!
  export DISPLAY="$XDISP"
  sleep 2
else
  XVFB_PID=""
fi

# --- monta o batch ---
{
  echo "genome hg38"
  shopt -s nullglob
  for t in "$EXT_DIR"/*.bw "$EXT_DIR"/*.bb "$EXT_DIR"/*.bed.gz; do
    [ -e "$t" ] && echo "load $(pwd)/$t"
  done
  echo "load $(pwd)/$BED_WITH"
  echo "load $(pwd)/$BED_WITHOUT"
  echo "snapshotDirectory $(pwd)/$OUTDIR"
} > igv_batch_option1.txt

awk -F'\t' '{print "goto "$1":"($2+1)"-"$3; print "snapshot"}' \
  "$BED_WITH" >> igv_batch_option1.txt
echo "exit" >> igv_batch_option1.txt

# --- roda IGV em batch (headless) ---
micromamba run -n igv igv -b "$(pwd)/igv_batch_option1.txt"

# --- encerra Xvfb somente se nós subimos ---
if [ -n "$XVFB_PID" ]; then kill "$XVFB_PID" 2>/dev/null || true; fi

echo "PNGs gerados em: $(pwd)/$OUTDIR"
