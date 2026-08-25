#!/usr/bin/env bash
# Gera snapshots IGV: ANOTACAO APENAS DA AMOSTRA COM VARIANTE + READS (BAMs) DE AMBAS
# as amostras (variante + controle) por locus.
#
# Requisitos:
#   - TSV gerado por str_samples_to_bed.R: str_samples_bams.tsv
#       col: gene  STRs_ID  chr  start0  end  variant_sample  variant_bam  control_sample  control_bam
#   - str_samples_with_variant.bed (anotacao da variante, 9 colunas)
#   - external_tracks/ (sinais opcionais .bw/.bb/.bed.gz)
#   - env micromamba 'igv' ativo (igv no PATH) ou igv no PATH
#
# Uso no cluster:
#   micromamba activate igv
#   Rscript str_samples_to_bed.R        # (re)gera str_samples_bams.tsv
#   bash run_igv_paired_reads.sh
#
# Saida: snapshots_paired/<gene>.png
set -u

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# ---------- resolver comando igv ----------
if command -v igv >/dev/null 2>&1; then
  IGV_CMD="igv"
elif command -v micromamba >/dev/null 2>&1; then
  IGV_CMD="$(micromamba run -n igv which igv 2>/dev/null)"
  [ -z "$IGV_CMD" ] && IGV_CMD="$(micromamba shell activate igv >/dev/null 2>&1; command -v igv)"
fi
if [ -z "${IGV_CMD:-}" ] || ! command -v "$IGV_CMD" >/dev/null 2>&1; then
  echo "ERRO: 'igv' nao encontrado. Ative o env 'igv' (micromamba activate igv)." >&2
  exit 1
fi
echo "Usando IGV: $IGV_CMD"

# ---------- garantir Xvfb :99 vivo ----------
XSOCK="/tmp/.X11-unix/X99"
ensure_xvfb() {
  if [ -S "$XSOCK" ] && command -v xdpyinfo >/dev/null 2>&1 && DISPLAY=:99 xdpyinfo >/dev/null 2>&1; then
    return 0
  fi
  command -v Xvfb >/dev/null 2>&1 || { echo "ERRO: Xvfb ausente." >&2; exit 1; }
  pkill -f "Xvfb :99" 2>/dev/null; sleep 1
  Xvfb :99 -screen 0 1280x1024x24 >/tmp/xvfb_igv.log 2>&1 &
  for i in $(seq 1 20); do
    [ -S "$XSOCK" ] && DISPLAY=:99 xdpyinfo >/dev/null 2>&1 && return 0
    sleep 0.5
  done
  echo "ERRO: falhou ao subir Xvfb :99." >&2; exit 1
}
ensure_xvfb
export DISPLAY=:99

# ---------- checagens ----------
TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)." >&2; exit 1; }
[ -f "$ANN" ] || { echo "ERRO: $ANN ausente." >&2; exit 1; }

mkdir -p snapshots_paired

# ---------- loop por locus ----------
first=1
while IFS=$'\t' read -r gene strs chr start0 end vsample vbam csample cbam; do
  # pula cabecalho
  if [ "$first" -eq 1 ]; then first=0; [ "$gene" = "gene" ] && continue; fi
  [ -z "$gene" ] && continue

  echo "==> $gene ($strs) chr$chr:start0=$start0 end=$end"
  echo "    variante: $vsample | controle: $csample"

  batch="igv_batch_${gene}.txt"
  {
    echo "genome hg38"
    # sinais opcionais
    for t in external_tracks/*.bw external_tracks/*.bb external_tracks/*.bed.gz; do
      [ -e "$t" ] && echo "load $(pwd)/$t"
    done
    # anotacao APENAS da amostra com variante
    echo "load $(pwd)/$ANN"
    # reads de ambas as amostras
    if [ -f "$vbam" ]; then echo "load $vbam"; else echo "WARN: BAM variante ausente: $vbam" >&2; fi
    if [ -f "$cbam" ]; then echo "load $cbam"; else echo "WARN: BAM controle ausente: $cbam" >&2; fi
    echo "snapshotDirectory $(pwd)/snapshots_paired"
  } > "$batch"

  # 1-based IGV
  s1=$(( start0 + 1 ))
  echo "goto ${chr}:${s1}-${end}" >> "$batch"
  echo "snapshot ${gene}.png" >> "$batch"
  echo "exit" >> "$batch"

  "$IGV_CMD" -b "$(pwd)/$batch"
  rm -f "$batch"
done < "$TSV"

echo "Pronto. Snapshots em: $(pwd)/snapshots_paired/"
