#!/usr/bin/env bash
# IGV.js (navegador) para a variante NKAIN2.
# Anotacao so da amostra com variante + reads (BAMs) das duas amostras.
# Os BAMs sao EXTRAIDOS para a regiao do locus (+/- FLANK). Tenta regiao
# (rapido, precisa de indice); se falhar, faz streaming + filtro por coordenada.
# O contig e renomeado p/ "chrN" (padrao hg38 do IGV.js) caso o BAM original
# use nomenclatura sem "chr" (ex.: "10" em vez de "chr10").
# Uso: bash igv_per_variant/NKAIN2.sh   (requer samtools + str_samples_bams.tsv)
set -u

GENE="NKAIN2"; PORT="8204"; FLANK=1000
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$BASE"

TSV="str_samples_bams.tsv"
ANN="str_samples_with_variant.bed"
[ -f "$TSV" ] || { echo "ERRO: $TSV ausente (rode Rscript str_samples_to_bed.R)."; exit 1; }
[ -f "$ANN" ] || echo "WARN: $ANN ausente; anotacao pode faltar." >&2

row="$(awk -F'\t' -v g="$GENE" '$1==g' "$TSV")"
[ -z "$row" ] && { echo "ERRO: $GENE nao encontrado no TSV."; exit 1; }
chr="$(echo "$row" | cut -f3)"; start0="$(echo "$row" | cut -f4)"; end="$(echo "$row" | cut -f5)"
vbam="$(echo "$row" | cut -f7 | sed 's|//*|/|g')"
cbam="$(echo "$row" | cut -f9 | sed 's|//*|/|g')"
s1=$((start0 + 1))

command -v samtools >/dev/null 2>&1 || { echo "ERRO: samtools ausente. Instale: micromamba install -n igv -c bioconda samtools"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "ERRO: python ausente."; exit 1; }

OUT="/tmp/igvjs_${GENE}"; mkdir -p "$OUT"

# BED so deste locus (fallback: BED inteiro)
awk -v c="$chr" -v s="$start0" -v e="$end" '($1==c && $2==s && $3==e)' "$ANN" > "$OUT/$GENE.bed"
if [ -s "$OUT/$GENE.bed" ]; then BED_URL="$OUT/$GENE.bed"; else BED_URL="$(pwd)/$ANN"; fi
BED_URL="$(echo "$BED_URL" | sed 's|//*|/|g')"

# contig como aparece no HEADER do BAM
bam_chr() {
  local bam="$1" want="$2" hit=""
  hit="$(samtools view -H "$bam" | awk -v w="$want" -F'\t' '{for(i=1;i<=NF;i++) if($i=="SN:"w){print w; exit}}')"
  if [ -n "$hit" ]; then echo "$want"; return; fi
  local alt; case "$want" in chr*) alt="${want#chr}";; *) alt="chr$want";; esac
  hit="$(samtools view -H "$bam" | awk -v w="$alt" -F'\t' '{for(i=1;i<=NF;i++) if($i=="SN:"w){print w; exit}}')"
  if [ -n "$hit" ]; then echo "$alt"; else echo "$want"; fi
}

# renomeia o contig do BAM extraido para $chr (padrao hg38), se necessario
norm_chr() {
  local bam="$1" cur="$2" want="$3"
  if [ "$cur" = "$want" ]; then
    samtools index "$bam" 2>>"/tmp/igvjs_${GENE}_samtools.log"
    return
  fi
  samtools view -h "$bam" \
    | awk -v cur="$cur" -v want="$want" 'BEGIN{FS=OFS="\t"} {
         if ($1 ~ /^@SQ/) sub("SN:"cur"\t","SN:"want"\t",$0);
         else if ($1 !~ /^@/ && $3==cur) $3=want;
         print }' \
    | samtools view -b - > "${bam}.tmp" 2>>"/tmp/igv_${GENE}_samtools.log" \
    && mv "${bam}.tmp" "$bam"
  samtools index "$bam" 2>>"/tmp/igv_${GENE}_samtools.log"
}

extract_bam() {
  local full="$1" label="$2"
  [ -f "$full" ] || { echo "WARN: BAM ausente: $full" >&2; return; }
  local cur; cur="$(bam_chr "$full" "$chr")"
  local fstart=$(( s1 - FLANK )); [ "$fstart" -lt 1 ] && fstart=1
  local fend=$(( end + FLANK ))
  local outb="$OUT/$GENE.$label.bam"
  echo "Extraindo $label (contig '$cur'): ${cur}:${fstart}-${fend}" >&2
  # 1) tenta regiao (rapido)
  if samtools view -b -h "$full" "${cur}:${fstart}-${fend}" > "$outb" 2>>"/tmp/igv_${GENE}_samtools.log"; then
    if [ -s "$outb" ]; then
      norm_chr "$outb" "$cur" "$chr"
      echo "$outb $outb.bai"; return
    fi
  fi
  # 2) fallback: streaming + filtro por coordenada (nao precisa de indice)
  echo "WARN: regiao falhou p/ $full; usando streaming (pode demorar)..." >&2
  samtools view -h "$full" \
    | awk -v c="$cur" -v a="$fstart" -v b="$fend" 'BEGIN{FS=OFS="\t"} {
         if ($1 ~ /^@/) { print; next }
         if ($3==c && $4>=a && $4<=b) print }' \
    | samtools view -b - > "$outb" 2>>"/tmp/igv_${GENE}_samtools.log"
  if [ -s "$outb" ]; then
    norm_chr "$outb" "$cur" "$chr"
    echo "$outb $outb.bai"
  else
    echo "WARN: extract vazio para $full" >&2
  fi
}
read -r vb_out vb_idx <<< "$(extract_bam "$vbam" variant)"
read -r cb_out cb_idx <<< "$(extract_bam "$cbam" control)"

# tracks.json (sem virgula inicial; paste insere o separador)
tt="$(mktemp)"
echo "{\"type\":\"annotation\",\"name\":\"variante $GENE\",\"url\":\"$BED_URL\",\"format\":\"bed\"}" > "$tt"
add_extracted() {
  local b="$1" idx="$2" name="$3"
  if [ -n "$b" ] && [ -f "$b" ] && [ -n "$idx" ] && [ -f "$idx" ]; then
    echo "{\"type\":\"alignment\",\"name\":\"$name\",\"url\":\"$b\",\"indexURL\":\"$idx\"}" >> "$tt"
  else
    echo "WARN: track nao adicionado: $b ($idx)" >&2
  fi
}
add_extracted "$vb_out" "$vb_idx" "$(basename "$vbam")"
add_extracted "$cb_out" "$cb_idx" "$(basename "$cbam")"
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

# Servidor HTTP a partir da raiz (/) com CORS. Sem kwargs directory/bind (compat. Python antigo).
python - "$PORT" <<'PY' >"/tmp/igvjs_${GENE}.log" 2>&1 &
import sys, http.server, socketserver, os
port = int(sys.argv[1])
os.chdir('/')
class H(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()
socketserver.TCPServer.allow_reuse_address = True
httpd = socketserver.TCPServer(('', port), H)
httpd.serve_forever()
PY
PID=$!
trap "kill $PID 2>/dev/null" EXIT

sleep 1
if command -v curl >/dev/null 2>&1; then
  echo "Auto-teste (HTTP 200 esperado):"
  for u in "$BED_URL" "$vb_out" "$vb_idx" "$cb_out" "$cb_idx"; do
    [ -n "$u" ] && code=$(curl -s -o /dev/null -w '%{http_code}' "http://localhost:$PORT$u") && echo "  [$code]  $u"
  done
else
  echo "WARN: curl ausente; auto-teste ignorado."
fi

echo
echo "============================================================"
echo "Variante $GENE"
echo "Regiao: ${chr}:${s1}-${end} (flank +/-${FLANK})"
echo "BAM variante extraido: $vb_out"
echo "BAM controle extraido: $cb_out"
echo "Abra no navegador:  http://localhost:$PORT/tmp/igvjs_$GENE/index.html"
echo "No PC:  ssh -L $PORT:localhost:$PORT Carlos_Chagas"
echo "Log IGV: /tmp/igvjs_$GENE.log   samtools: /tmp/igv_${GENE}_samtools.log"
echo "============================================================"

wait "$PID"
