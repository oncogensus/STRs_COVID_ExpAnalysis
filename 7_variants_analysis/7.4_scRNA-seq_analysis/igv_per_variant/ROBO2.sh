#!/usr/bin/env bash
# IGV.js (navegador) para a variante ROBO2.
# Anotacao so da amostra com variante + BAMs das duas amostras.
# Uso: bash igv_per_variant/ROBO2.sh   (requer str_samples_bams.tsv + str_samples_with_variant.bed)
set -u

GENE="ROBO2"; PORT="8201"
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)."; exit 1; }
[ -f "$ANN" ] || echo "WARN: $ANN ausente; anotacao pode faltar." >&2

row="$(awk -F'\t' -v g="$GENE" '$1==g' "$TSV")"
[ -z "$row" ] && { echo "ERRO: $GENE nao encontrado no TSV."; exit 1; }
chr="$(echo "$row" | cut -f3)"; start0="$(echo "$row" | cut -f4)"; end="$(echo "$row" | cut -f5)"
vbam="$(echo "$row" | cut -f7)"; cbam="$(echo "$row" | cut -f9)"
s1=$((start0 + 1))

DATA="$BASE/igv_per_variant/data/$GENE"; mkdir -p "$DATA"

# BED so deste locus (fallback: BED inteiro)
awk -v c="$chr" -v s="$start0" -v e="$end" '($1==c && $2==s && $3==e)' "$ANN" > "$DATA/$GENE.bed"
if [ ! -s "$DATA/$GENE.bed" ]; then ln -sf "$(pwd)/$ANN" "$DATA/$GENE.bed"; fi

# symlinks dos BAMs (+ .bai)
for b in "$vbam" "$cbam"; do
  if [ -f "$b" ]; then
    ln -sf "$b" "$DATA/$(basename "$b")"
    [ -f "$b.bai" ] && ln -sf "$b.bai" "$DATA/$(basename "$b").bai"
  else
    echo "WARN: BAM ausente: $b" >&2
  fi
done

# tracks.json
tt="$(mktemp)"
echo "{\"type\":\"annotation\",\"name\":\"variante $GENE\",\"url\":\"$GENE.bed\",\"format\":\"bed\"}" > "$tt"
[ -f "$vbam" ] && echo "{\"type\":\"alignment\",\"name\":\"$(basename "$vbam")\",\"url\":\"$(basename "$vbam")\",\"indexURL\":\"$(basename "$vbam").bai\"}" >> "$tt"
[ -f "$cbam" ] && echo "{\"type\":\"alignment\",\"name\":\"$(basename "$cbam")\",\"url\":\"$(basename "$cbam")\",\"indexURL\":\"$(basename "$cbam").bai\"}" >> "$tt"
{ echo "["; paste -sd, "$tt"; echo "]"; } > "$DATA/tracks.json"
rm -f "$tt"

# index.html (igv.js via CDN)
sed "s/__LOCUS__/${chr}:${s1}-${end}/" > "$DATA/index.html" <<'HTML'
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
python -m http.server "$PORT" --directory "$DATA" >"/tmp/igvjs_${GENE}.log" 2>&1 &
PID=$!
trap "kill $PID 2>/dev/null" EXIT

echo
echo "============================================================"
echo "Variante $GENE  ->  http://localhost:$PORT"
echo "No PC:  ssh -L $PORT:localhost:$PORT Carlos_Chagas"
echo "Abrir no navegador: http://localhost:$PORT"
echo "Regiao: ${chr}:${s1}-${end}"
echo "============================================================"

wait "$PID"
