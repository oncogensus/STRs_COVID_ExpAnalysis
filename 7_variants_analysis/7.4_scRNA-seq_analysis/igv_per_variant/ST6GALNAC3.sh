#!/usr/bin/env bash
# IGV.js (navegador) para a variante ROBO2.
# Anotacao so da amostra com variante + BAMs das duas amostras.
# Abordagem: o servidor HTTP sobe a PARTIR DA RAIZ (/) e os tracks usam
# caminhos absolutos, evitando problemas de symlink fora do diretorio.
# Uso: bash igv_per_variant/ROBO2.sh   (requer str_samples_bams.tsv + str_samples_with_variant.bed)
set -u

GENE="ST6GALNAC3"; PORT="8208"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)."; exit 1; }
[ -f "$ANN" ] || echo "WARN: $ANN ausente; anotacao pode faltar." >&2

row="$(awk -F'\t' -v g="$GENE" '$1==g' "$TSV")"
[ -z "$row" ] && { echo "ERRO: $GENE nao encontrado no TSV."; exit 1; }
chr="$(echo "$row" | cut -f3)"; start0="$(echo "$row" | cut -f4)"; end="$(echo "$row" | cut -f5)"
# normaliza barras duplas (ex.: recal//x.bam)
vbam="$(echo "$row" | cut -f7 | sed 's|//*|/|g')"
cbam="$(echo "$row" | cut -f9 | sed 's|//*|/|g')"
s1=$((start0 + 1))

OUT="/tmp/igvjs_${GENE}"; mkdir -p "$OUT"

# BED so deste locus (fallback: BED inteiro, por caminho absoluto)
awk -v c="$chr" -v s="$start0" -v e="$end" '($1==c && $2==s && $3==e)' "$ANN" > "$OUT/$GENE.bed"
if [ -s "$OUT/$GENE.bed" ]; then BED_URL="$OUT/$GENE.bed"; else BED_URL="$(pwd)/$ANN"; fi
# normaliza barras duplas em todas as URLs
BED_URL="$(echo "$BED_URL" | sed 's|//*|/|g')"
vbam="$(echo "$vbam" | sed 's|//*|/|g')"
cbam="$(echo "$cbam" | sed 's|//*|/|g')"

# helper: descobre o arquivo de indice (.bam.bai ou .bai)
idx_of() {
  local b="$1" c=""
  for cand in "$b.bai" "${b%.bam}.bai"; do
    [ -f "$cand" ] && { c="$cand"; break; }
  done
  [ -n "$c" ] && printf '%s' "$c"
}

tt="$(mktemp)"
echo "{\"type\":\"annotation\",\"name\":\"variante $GENE\",\"url\":\"$BED_URL\",\"format\":\"bed\"}" > "$tt"

add_bam() {
  local b="$1"
  if [ -f "$b" ]; then
    local iu; iu="$(idx_of "$b")"
    if [ -n "$iu" ]; then
      echo "{\"type\":\"alignment\",\"name\":\"$(basename "$b")\",\"url\":\"$b\",\"indexURL\":\"$iu\"}" >> "$tt"
    else
      echo "WARN: indice .bai ausente para $b" >&2
    fi
  else
    echo "WARN: BAM ausente: $b" >&2
  fi
}
add_bam "$vbam"
add_bam "$cbam"

{ echo "["; paste -sd, "$tt"; echo "]"; } > "$OUT/tracks.json"
rm -f "$tt"

# index.html (igv.js via CDN)
sed "s/__LOCUS__/${chr}:${s1}-${end}/" > "$OUT/index.html" <<'HTML'
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

command -v python >/dev/null 2>&1 || { echo "ERRO: python ausente."; exit 1; }
python -m http.server "$PORT" --directory / >"/tmp/igvjs_${GENE}.log" 2>&1 &
PID=$!
trap "kill $PID 2>/dev/null" EXIT

echo
echo "============================================================"
echo "Variante $GENE"
echo "Abra no navegador:  http://localhost:$PORT/tmp/igvjs_$GENE/index.html"
echo "No PC:  ssh -L $PORT:localhost:$PORT Carlos_Chagas"
echo "Regiao: ${chr}:${s1}-${end}"
echo "============================================================"

wait "$PID"
