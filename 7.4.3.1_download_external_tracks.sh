#!/bin/bash
# ==============================================================================
# 7.4.3.1_download_external_tracks.sh
# ==============================================================================
# Downloads external genomic tracks required by 7.4.3.2_trackviewer_STR_viz.ipynb
#
# Output: all files saved under ./external_tracks/
# ==============================================================================
set -euo pipefail

OUTDIR="external_tracks"
mkdir -p "$OUTDIR"

LOGFILE="$OUTDIR/download_external_tracks.log"
exec > >(
    while IFS= read -r line; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') $line"
    done | tee -a "$LOGFILE"
) 2>&1

echo "=============================================="
echo " download started — log: $LOGFILE"
echo "=============================================="

cd "$OUTDIR"

# TFs that failed to download in B16 (best-effort); reported in the final summary.
declare -a FAILED_TFS=()

# ------------------------------------------------------------------------------
# Helper 1: standard HTTP download with HTML error validation
# ------------------------------------------------------------------------------
download() {
  local url="$1" out="$2"; shift 2
  local fallbacks=("$@")
  if [ -f "$out" ]; then
    echo "   [skip] $out already present"
    return 0
  fi
  local all=("$url" "${fallbacks[@]}")
  local u http_code start_ts end_ts size
  for u in "${all[@]}"; do
    echo "   GET  $u"
    start_ts=$(date +%s)
    http_code=$(curl -sSL -o "$out.tmp" -w '%{http_code}' --retry 2 --retry-delay 3 "$u") || http_code="000"
    end_ts=$(date +%s)
    if [ "$http_code" = "200" ] && [ -s "$out.tmp" ]; then
      # magic-byte sanity check: reject HTML/error pages served with a 200
      if head -c 200 "$out.tmp" | grep -qiE '<html|<!DOCTYPE'; then
        echo "   [FAIL] HTTP $http_code but payload is an HTML/error page; removed"
        rm -f "$out.tmp"
        continue
      fi
      size=$(stat -c%s "$out.tmp" 2>/dev/null || wc -c < "$out.tmp")
      mv "$out.tmp" "$out"
      echo "   [OK]   HTTP $http_code — ${size} bytes, $((end_ts - start_ts))s -> $out"
      return 0
    fi
    echo "   [FAIL] HTTP $http_code ($((end_ts - start_ts))s) for: $u"
    rm -f "$out.tmp"
  done
  echo "   [FAIL] all $(( ${#all[@]} )) candidate URL(s) failed for: $out"
  return 1
}

# ------------------------------------------------------------------------------
# Helper 2: ReMap URL candidates lookup
# ------------------------------------------------------------------------------
remap_tf_candidates() {
  local tf="$1"
  local lc; lc=$(printf '%s' "$tf" | tr '[:upper:]' '[:lower:]')
  local tc; tc=$(printf '%s' "$tf" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')
  
  echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$tf/remap2022_${tf}_all_macs2_hg38_v1_0.bed.gz"
  [ "$lc" != "$tf" ] && echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$lc/remap2022_${lc}_all_macs2_hg38_v1_0.bed.gz"
  [ "$tc" != "$tf" ] && echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$tc/remap2022_${tc}_all_macs2_hg38_v1_0.bed.gz"
}

# ------------------------------------------------------------------------------
# Helper 3: CLI Fallbacks (1. UCSC MySQL Public | 2. JASPAR 2026 bigBed)
# ------------------------------------------------------------------------------
extract_tf_cli_fallbacks() {
  local tf="$1"
  local out="$2"

  export -f extract_tf_cli_fallbacks 2>/dev/null || true

  if [ -f "$out" ]; then
    return 0
  fi

  echo "   [fallback] ReMap failed for $tf. Trying CLI Method 1: UCSC MySQL..."
  if command -v mysql &>/dev/null; then
    mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -D hg38 \
      -e "SELECT chrom, chromStart, chromEnd, name, score FROM encRegTfbsClustered WHERE name LIKE '%${tf}%';" 2>/dev/null \
      | sed '1d' | gzip -c > "$out.tmp" || true
    
    if [ -s "$out.tmp" ]; then
      mv "$out.tmp" "$out"
      echo "   [OK]   Successfully extracted $tf via UCSC MySQL -> $out"
      return 0
    fi
    rm -f "$out.tmp"
  else
    echo "   [info] mysql client not found in PATH; skipping UCSC MySQL extraction."
  fi

  echo "   [fallback] Trying CLI Method 2: Extracting $tf from JASPAR2026.bb..."
  if [ -f "JASPAR2026.bb" ]; then
    # Ensure bigBedToBed executable is available
    local b2b_cmd="bigBedToBed"
    if ! command -v bigBedToBed &>/dev/null; then
      if [ ! -x "./bigBedToBed" ]; then
        echo "   [info] Fetching standalone bigBedToBed binary from UCSC..."
        curl -sSL -o ./bigBedToBed "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed" 2>/dev/null || true
        chmod +x ./bigBedToBed 2>/dev/null || true
      fi
      b2b_cmd="./bigBedToBed"
    fi

    if [ -x "$(command -v "$b2b_cmd")" ] || [ -x "$b2b_cmd" ]; then
      "$b2b_cmd" JASPAR2026.bb stdout 2>/dev/null | grep -i "${tf}" | gzip -c > "$out.tmp" || true
      if [ -s "$out.tmp" ]; then
        mv "$out.tmp" "$out"
        echo "   [OK]   Successfully extracted $tf via JASPAR2026.bb -> $out"
        return 0
      fi
      rm -f "$out.tmp"
    else
      echo "   [FAIL] Could not run bigBedToBed."
    fi
  else
    echo "   [FAIL] JASPAR2026.bb file not present for local extraction."
  fi

  return 1
}

echo "=============================================="
echo " Downloading external tracks -> $OUTDIR"
echo "=============================================="

# B1. RepeatMasker
echo "[B1] RepeatMasker hg38 (UCSC rmsk)"
download "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" "rmsk.hg38.bed.gz"

# B2. GTEx eQTL
echo "[B2] GTEx eQTL SNP list"
if [ ! -f gtex_eqTL_rsids.tsv ]; then
  cat > gtex_eqTL_rsids.tsv <<'EOF'
rsID	tissue
rs17016404	Brain_Hypothalamus
rs17016416	Brain_Hypothalamus
rs78303965	Brain_Hypothalamus
rs75109825	Brain_Hypothalamus
rs17016490	Brain_Hypothalamus
rs11167599	Brain_Hypothalamus
rs2968279	Brain_Cerebellum
rs2920747	Brain_Cerebellum
rs72890312	Cells_Cultured_fibroblasts
rs73114451	Cells_Cultured_fibroblasts
rs7311445	Cells_Cultured_fibroblasts
rs12526162	Brain_Cortex
rs74471882	Brain_Cortex
rs4708037	Brain_Cerebellar_Hemisphere
rs34353921	Brain_Cerebellar_Hemisphere
rs71574845	Brain_Frontal_Cortex_BA9
rs2753166	Brain_Frontal_Cortex_BA9
rs13220639	Brain_Frontal_Cortex_BA9
rs138555200	Brain_Hippocampus
EOF
fi
echo "   -> saved gtex_eqTL_rsids.tsv ($(wc -l < gtex_eqTL_rsids.tsv) lines)"

# B3. CTCF
echo "[B3] CTCF ChIP-seq (ENCODE brain, ENCFF910VLV)"
download "https://www.encodeproject.org/files/ENCFF910VLV/@@download/ENCFF910VLV.bigWig" "CTCF_ENCFF910VLV.bw"

# B4. DNase-seq
echo "[B4] DNase-seq accessibility (ENCODE brain, ENCFF110ZAI)"
download "https://www.encodeproject.org/files/ENCFF110ZAI/@@download/ENCFF110ZAI.bigWig" "DNase_brain.bw"

# B5. H3K27ac
echo "[B5] H3K27ac (ENCODE brain, ENCFF542HTC)"
download "https://www.encodeproject.org/files/ENCFF542HTC/@@download/ENCFF542HTC.bigWig" "H3K27ac_brain.bw"

# B6. H3K27me3
echo "[B6] H3K27me3 (ENCODE brain, ENCFF722UXM)"
download "https://www.encodeproject.org/files/ENCFF722UXM/@@download/ENCFF722UXM.bigWig" "H3K27me3_brain.bw"

# B7. Main ReMap TFs
echo "[B7] ReMap TF ChIP-seq peaks"
declare -A TF_FILES=(
  [CTCF]="remap2022_ctcf_all_macs2_hg38_v1_0.bed.gz"
  [REST]="remap2022_rest_all_macs2_hg38_v1_0.bed.gz"
  [KLF9]="remap2022_klf9_all_macs2_hg38_v1_0.bed.gz"
  [GATA2]="remap2022_gata2_all_macs2_hg38_v1_0.bed.gz"
  [ZNF384]="remap2022_znf384_all_macs2_hg38_v1_0.bed.gz"
  [MNT]="remap2022_mnt_all_macs2_hg38_v1_0.bed.gz"
)
for tf in "${!TF_FILES[@]}"; do
  echo "   $tf"
  mapfile -t cands < <(remap_tf_candidates "$tf")
  download "${cands[0]}" "${TF_FILES[$tf]}" "${cands[@]:1}"
done

# B8. ReMap Density
echo "[B8] ReMap 2022 density (hg38)"
download "https://remap.univ-amu.fr/storage/public/hubReMap2022/hg38/bigBed/hg38.bw" "remap2022_density_hg38.bw"

# B9. ClinGen
echo "[B9] ClinGen dosage sensitivity (GRCh38)"
download "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv" "ClinGen_gene_curation_list_GRCh38.tsv"

# B10. JARVIS
if [ "${DOWNLOAD_LARGE:-0}" = "1" ]; then
  echo "[B10] JARVIS non-coding constraint score (hg38)"
  download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/jarvis/jarvis.bw" "jarvis.bw"
else
  echo "[B10] SKIP JARVIS bigWig (21 GB). Notebook queries remote URL by window."
fi

# B11. UKB Depletion
echo "[B11] UK Biobank / deCODE Depletion Rank Score (hg38)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/ukbDepletion/ukbDepletion.bw" "ukbDepletion.bw"

# B12. ReMap Atlas UCSC
echo "[B12] ReMap Atlas density (hg38, UCSC)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMapDensity2022.bw" "remap_density_ucsc.bw"

# B13. JASPAR 2026 (Baixado antes do B16 para permitir extração local se necessário)
echo "[B13] JASPAR CORE TFBS (hg38, JASPAR2026.bb)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2026.bb" "JASPAR2026.bb"

# B14. DECIPHER CNVs
echo "[B14] DECIPHER common CNVs (hg38)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/dbVar/common_decipher.bb" "decipher_common.bb"

# B15. TRExplorer
echo "[B15] TRExplorer tandem repeat catalog (hg38)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/strVar/trexplorer.bb" "trexplorer.bb"

# B16. Extra ReMap TFs com Fallbacks CLI Automáticos
echo "[B16] ReMap TF ChIP-seq peaks (extra TFs)"
declare -A EXTRA_TF_FILES=(
  [NACC2]="remap2022_nacc2_all_macs2_hg38_v1_0.bed.gz"
  [FOXB1]="remap2022_foxb1_all_macs2_hg38_v1_0.bed.gz"
  [FOXK1]="remap2022_foxk1_all_macs2_hg38_v1_0.bed.gz"
  [PATZ1]="remap2022_patz1_all_macs2_hg38_v1_0.bed.gz"
  [RAD21]="remap2022_rad21_all_macs2_hg38_v1_0.bed.gz"
  [TCF4]="remap2022_tcf4_all_macs2_hg38_v1_0.bed.gz"
  [MITF]="remap2022_mitf_all_macs2_hg38_v1_0.bed.gz"
  [USF2]="remap2022_usf2_all_macs2_hg38_v1_0.bed.gz"
  [RAD51]="remap2022_rad51_all_macs2_hg38_v1_0.bed.gz"
  [GATA1]="remap2022_gata1_all_macs2_hg38_v1_0.bed.gz"
)

declare -a FAILED_TFS=()
for tf in "${!EXTRA_TF_FILES[@]}"; do
  echo "   $tf"
  mapfile -t cands < <(remap_tf_candidates "$tf")
  
  if download "${cands[0]}" "${EXTRA_TF_FILES[$tf]}" "${cands[@]:1}"; then
    :
  elif extract_tf_cli_fallbacks "$tf" "${EXTRA_TF_FILES[$tf]}"; then
    :
  else
    echo "   [WARN] $tf unavailable on ReMap 2022 and CLI fallbacks failed."
    FAILED_TFS+=("$tf")
  fi
done

# B17. ENCODE CTCF
echo "[B17] ENCODE wgEncodeReg4TfChip CTCF ChIP-seq (ENCFF341CQE)"
download "https://www.encodeproject.org/files/ENCFF341CQE/@@download/ENCFF341CQE.bigWig" "wgEncodeReg4TfChip_ENCFF341CQE.bw"

echo ""
echo "=============================================="
echo " DONE. Tracks available in: $(pwd)"
if [ ${#FAILED_TFS[@]} -gt 0 ]; then
  echo ""
  echo " [!] The following B16 TFs FAILED to download and are NOT present:"
  for t in "${FAILED_TFS[@]}"; do echo "     - $t"; done
  echo "     (NACC2 and FOXB1 are known to be absent from ReMap 2022 hg38;"
  echo "      NACC2 is queried via JASPAR2026.bb instead.)"
fi
ls -lh
echo "=============================================="