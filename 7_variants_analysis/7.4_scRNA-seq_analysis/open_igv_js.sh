#!/usr/bin/env bash
# open_igv_js.sh
# Abre as amostras (BED da variante + BAMs das duas amostras) no NAVEGADOR via IGV.js.
# Nao precisa de VNC. Requer apenas python (qualquer 3.7+) e internet no PC para o CDN do igv.js.
# O servidor HTTP sobe no cluster e e acessado por tunel SSH.
#
# Uso no cluster:
#   micromamba activate igv
#   Rscript str_samples_to_bed.R
#   bash open_igv_js.sh
# No PC:
#   ssh -L 8000:localhost:8000 Carlos_Chagas
#   abrir http://localhost:8000
set -u

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
DATA="$DIR/igvjs_data"
PORT="8000"

[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)." >&2; exit 1; }
[ -f "$ANN" ] || { echo "ERRO: $ANN ausente." >&2; exit 1; }

mkdir -p "$DATA"
rm -f "$DATA"/*.bam "$DATA"/*.bai "$DATA"/*.bed 2>/dev/null

# symlink do BED da variante
ln -sf "$(pwd)/$ANN" "$DATA/$ANN"

# symlink dos BAMs (+ .bai) e descoberta do primeiro locus
first_region=""
tracks_tmp="$(mktemp)"
echo "{\"type\":\"annotation\",\"name\":\"variante ($ANN)\",\"url\":\"$ANN\",\"format\":\"bed\"}" > "$tracks_tmp"

while IFS=$'\t' read -r gene strs chr start0 end vsample vbam csample cbam; do
  [ "$gene" = "gene" ] && continue
  [ -z "$gene" ] && continue
  for b in "$vbam" "$cbam"; do
    if [ -f "$b" ]; then
      bn="$(basename "$b")"
      ln -sf "$b" "$DATA/$bn"
      [ -f "$b.bai" ] && ln -sf "$b.bai" "$DATA/$bn.bai"
      echo "{\"type\":\"alignment\",\"name\":\"$bn\",\"url\":\"$bn\",\"indexURL\":\"$bn.bai\"}" >> "$tracks_tmp"
    else
      echo "WARN: BAM ausente: $b" >&2
    fi
  done
  if [ -z "$first_region" ]; then
    s1=$(( start0 + 1 )); first_region="${chr}:${s1}-${end}"
  fi
done < "$TSV"

# monta tracks.json (array JSON valido)
{ echo "["; paste -sd, "$tracks_tmp"; echo "]"; } > "$DATA/tracks.json"
rm -f "$tracks_tmp"

# index.html (igv.js via CDN)
sed "s/__LOCUS__/${first_region:-chr1:1-1000}/" > "$DATA/index.html" <<'HTML'
<!DOCTYPE html>
<html lang="pt-br">
<head>
<meta charset="utf-8"/>
<title>IGV STR</title>
<script src="https://cdn.jsdelivr.net/npm/igv@3/dist/igv.min.js"></script>
<style>html,body{margin:0;height:100%}#igv{width:100%;height:100%}</style>
</head>
<body>
<div id="igv"></div>
<script>
const locus = "__LOCUS__";
fetch('tracks.json')
  .then(r => r.json())
  .then(tracks => igv.createBrowser(document.getElementById('igv'),
        { genome: 'hg38', locus: locus, tracks: tracks }))
  .catch(e => { document.body.insertAdjacentHTML('beforeend',
        '<pre style="color:red">Erro: '+e+'</pre>'); });
</script>
</body>
</html>
HTML

# sobe servidor HTTP (Python >=3.7 suporta Range, necessario p/ BAM)
command -v python >/dev/null 2>&1 || { echo "ERRO: python ausente." >&2; exit 1; }
python -m http.server "$PORT" --directory "$DATA" >/tmp/igvjs.log 2>&1 &
SERVER_PID=$!
trap "kill $SERVER_PID 2>/dev/null" EXIT

echo
echo "============================================================"
echo "Servidor IGV.js no cluster (porta $PORT)."
echo "No seu PC:"
echo "  ssh -L $PORT:localhost:$PORT Carlos_Chagas"
echo "  abrir http://localhost:$PORT"
echo "  (primeiro locus: ${first_region:-<none>})"
echo "============================================================"

wait "$SERVER_PID"
