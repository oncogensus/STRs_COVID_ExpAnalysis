#!/usr/bin/env bash
# open_igv_web.sh
# Abre o IGV de forma interativa ACESSIVEL PELO NAVEGADOR (noVNC).
# Pipeline:
#   Xvfb :99  ->  x0vncserver (tigervnc) na porta 5901  ->  websockify + noVNC na 6080
# Dependencias (conda-forge, env 'igv'):
#   tigervnc   (fornece x0vncserver e Xvfb)
#   websockify
#   git + internet (para clonar o noVNC uma unica vez)
#
# Uso no cluster:
#   micromamba activate igv
#   Rscript str_samples_to_bed.R
#   bash open_igv_web.sh
# No PC:
#   ssh -L 6080:localhost:6080 Carlos_Chagas
#   abrir http://localhost:6080/vnc.html
set -u

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

VNC_DISPLAY=":99"
VNC_PORT="5901"
WEB_PORT="6080"
NOVNC_DIR="$DIR/noVNC"

# ---------- resolver comando igv ----------
if command -v igv >/dev/null 2>&1; then
  IGV_CMD="igv"
elif command -v micromamba >/dev/null 2>&1; then
  IGV_CMD="$(micromamba run -n igv which igv 2>/dev/null)"
  [ -z "$IGV_CMD" ] && IGV_CMD="$(micromamba shell activate igv >/dev/null 2>&1; command -v igv)"
fi
if [ -z "${IGV_CMD:-}" ] || ! command -v "$IGV_CMD" >/dev/null 2>&1; then
  echo "ERRO: 'igv' nao encontrado. Ative o env 'igv'." >&2; exit 1
fi

# ---------- Xvfb ----------
ensure_xvfb() {
  XSOCK="/tmp/.X11-unix/X${VNC_DISPLAY#:}"
  if [ -S "$XSOCK" ] && command -v xdpyinfo >/dev/null 2>&1 && DISPLAY="$VNC_DISPLAY" xdpyinfo >/dev/null 2>&1; then
    return 0
  fi
  command -v Xvfb >/dev/null 2>&1 || { echo "ERRO: Xvfb ausente (instale tigervnc)." >&2; exit 1; }
  pkill -f "Xvfb $VNC_DISPLAY" 2>/dev/null; sleep 1
  Xvfb "$VNC_DISPLAY" -screen 0 1280x1024x24 >/tmp/xvfb_igv.log 2>&1 &
  for i in $(seq 1 20); do [ -S "$XSOCK" ] && break; sleep 0.5; done
}
ensure_xvfb
export DISPLAY="$VNC_DISPLAY"

# ---------- VNC server (x0vncserver do tigervnc, ou x11vnc) ----------
cleanup() {
  pkill -f "x0vncserver" 2>/dev/null
  pkill -f "x11vnc" 2>/dev/null
  pkill -f "websockify" 2>/dev/null
}
trap cleanup EXIT

if command -v x0vncserver >/dev/null 2>&1; then
  x0vncserver -display "$VNC_DISPLAY" -rfbport "$VNC_PORT" -SecurityTypes None >/tmp/x0vnc.log 2>&1 &
  echo "x0vncserver (tigervnc) na porta $VNC_PORT."
elif command -v x11vnc >/dev/null 2>&1; then
  x11vnc -display "$VNC_DISPLAY" -bg -nopw -listen localhost -rfbport "$VNC_PORT" >/tmp/x11vnc.log 2>&1 &
  echo "x11vnc na porta $VNC_PORT."
else
  echo "ERRO: nem x0vncserver nem x11vnc encontrados. Instale: micromamba install -n igv -c conda-forge tigervnc" >&2
  exit 1
fi
sleep 2
if ! command -v websockify >/dev/null 2>&1; then
  echo "ERRO: websockify ausente. Instale: micromamba install -n igv -c conda-forge websockify" >&2
  exit 1
fi

# ---------- noVNC (cliente web) ----------
if [ ! -d "$NOVNC_DIR" ]; then
  echo "Clonando noVNC (uma vez)..."
  git clone --depth 1 https://github.com/novnc/noVNC.git "$NOVNC_DIR" || {
    echo "ERRO: falhou clonar noVNC (precisa de git+internet)." >&2; exit 1; }
fi
websockify -D --web "$NOVNC_DIR" "$WEB_PORT" "localhost:$VNC_PORT"
echo "websockify na porta $WEB_PORT (web em $NOVNC_DIR)."

# ---------- checagens dos dados ----------
TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)." >&2; exit 1; }
[ -f "$ANN" ] || { echo "ERRO: $ANN ausente." >&2; exit 1; }

# ---------- batch IGV ----------
batch="igv_open_web.txt"
: > "$batch"
{
  echo "genome hg38"
  for t in external_tracks/*.bw external_tracks/*.bb external_tracks/*.bed.gz; do
    [ -e "$t" ] && echo "load $(pwd)/$t"
  done
  echo "load $(pwd)/$ANN"
} >> "$batch"

first_region=""
while IFS=$'\t' read -r gene strs chr start0 end vsample vbam csample cbam; do
  [ "$gene" = "gene" ] && continue
  [ -z "$gene" ] && continue
  [ -f "$vbam" ] && echo "load $vbam" || echo "WARN: BAM variante ausente: $vbam" >&2
  [ -f "$cbam" ] && echo "load $cbam" || echo "WARN: BAM controle ausente: $cbam" >&2
  if [ -z "$first_region" ]; then
    s1=$(( start0 + 1 )); first_region="${chr}:${s1}-${end}"
  fi
done < "$TSV" >> "$batch"

[ -n "$first_region" ] && echo "goto $first_region" >> "$batch"

echo
echo "============================================================"
echo "IGV aberto no display virtual. Acesse pelo NAVEGADOR:"
echo "  1) No seu PC:  ssh -L $WEB_PORT:localhost:$WEB_PORT Carlos_Chagas"
echo "  2) Navegador:  http://localhost:$WEB_PORT/vnc.html"
echo "  (primeiro locus: ${first_region:-<none>})"
echo "============================================================"

"$IGV_CMD" -b "$(pwd)/$batch"
rm -f "$batch"
