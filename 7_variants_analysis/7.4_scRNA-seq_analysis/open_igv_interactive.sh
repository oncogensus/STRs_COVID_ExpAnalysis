#!/usr/bin/env bash
# open_igv_interactive.sh
# Abre o IGV (modo interativo / janela) com:
#   - sinais de external_tracks/ (opcional)
#   - str_samples_with_variant.bed  (anotacao da variante)
#   - BAMs das DUAS amostras (variante + controle) de TODOS os loci
# e posiciona no primeiro locus. NAO gera snapshots nem sai.
#
# Para ver a janela no cluster, conecte com X forwarding:  ssh -Y usuario@cluster
# (ness caso o DISPLAY ja vem ajustado e o Xvfb nao sera usado).
# Se nao houver DISPLAY, sobe Xvfb :99 (display virtual; so util se houver VNC).
#
# Uso:
#   micromamba activate igv
#   bash open_igv_interactive.sh
set -u

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# carrega tambem o BED da amostra SEM variante? (0=nao, 1=sim)
LOAD_WITHOUT_BED=0

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

# ---------- DISPLAY ----------
if [ -z "${DISPLAY:-}" ]; then
  echo "DISPLAY vazio: subindo Xvfb :99 (display virtual)."
  XSOCK="/tmp/.X11-unix/X99"
  command -v Xvfb >/dev/null 2>&1 || { echo "ERRO: Xvfb ausente." >&2; exit 1; }
  pkill -f "Xvfb :99" 2>/dev/null; sleep 1
  Xvfb :99 -screen 0 1280x1024x24 >/tmp/xvfb_igv.log 2>&1 &
  for i in $(seq 1 20); do
    [ -S "$XSOCK" ] && break; sleep 0.5
  done
  export DISPLAY=:99
else
  echo "Usando DISPLAY=$DISPLAY (X forwarding?)."
fi

# ---------- checagens ----------
TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)." >&2; exit 1; }
[ -f "$ANN" ] || { echo "ERRO: $ANN ausente." >&2; exit 1; }

batch="igv_open_interactive.txt"
: > "$batch"
{
  echo "genome hg38"
  for t in external_tracks/*.bw external_tracks/*.bb external_tracks/*.bed.gz; do
    [ -e "$t" ] && echo "load $(pwd)/$t"
  done
  echo "load $(pwd)/$ANN"
  if [ "$LOAD_WITHOUT_BED" -eq 1 ] && [ -f str_samples_without_variant.bed ]; then
    echo "load $(pwd)/str_samples_without_variant.bed"
  fi
} >> "$batch"

# BAMs de todos os loci (variante + controle)
first_region=""
while IFS=$'\t' read -r gene strs chr start0 end vsample vbam csample cbam; do
  [ "$gene" = "gene" ] && continue
  [ -z "$gene" ] && continue
  if [ -f "$vbam" ]; then echo "load $vbam"; else echo "WARN: BAM variante ausente: $vbam" >&2; fi
  if [ -f "$cbam" ]; then echo "load $cbam"; else echo "WARN: BAM controle ausente: $cbam" >&2; fi
  if [ -z "$first_region" ]; then
    s1=$(( start0 + 1 ))
    first_region="${chr}:${s1}-${end}"
  fi
done < "$TSV" >> "$batch"

if [ -n "$first_region" ]; then
  echo "goto $first_region" >> "$batch"
fi
# SEM 'exit' -> IGV abre a janela e fica interativo

echo "Carregando IGV (janela interativa). Primeiro locus: ${first_region:-<none>}"
"$IGV_CMD" -b "$(pwd)/$batch"
rm -f "$batch"
