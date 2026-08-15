#!/bin/bash
# ==============================================================================
# 7.4.3.1_download_external_tracks.sh
# ==============================================================================
# Downloads the external genomic tracks required by 7.4.3.2_trackviewer_STR_viz.ipynb
# to visualize the 8 filtered STRs from Relatorio_STR_Final_Integral.pdf
# with the trackViewer R package.
#
# Tracks downloaded (hg38 / GRCh38):
#   B1. RepeatMasker           (UCSC)
#   B2. GTEx eQTL SNPs         (GTEx Portal - selected significant variants)
#   B3. CTCF ChIP-seq (ENCODE) (brain, dorsolateral prefrontal cortex)
#   B4. DNase-seq              (ENCODE, brain)
#   B5. H3K27ac                (ENCODE, brain - dorsolateral prefrontal cortex)
#   B6. H3K27me3 (ENCODE)      (brain - repressive mark, replaces retired Roadmap)
#   B7. TF ChIP-seq (ReMap)    (CTCF, REST, KLF9, GATA2, ZNF384, MNT)
#   B8. ReMap Density          (ReMap 2022 track hub)
#   B9. ClinGen dosage sensitivity (gene curation list, GRCh38)
#
# Usage:
#   bash 7.4.3.1_download_external_tracks.sh
#
# Output: all files saved under ./external_tracks/
#
# NOTE: every download is validated after retrieval; HTML error pages (e.g. a
# 404 response) are detected and removed so the notebook never imports a bogus
# file. Only tracks that were actually downloaded will be shown.
# ==============================================================================
set -euo pipefail

OUTDIR="external_tracks"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# TFs that failed to download in B16 (best-effort); reported in the final summary.
declare -a FAILED_TFS=()

# ------------------------------------------------------------------------------
# Helper: fetch a file, trying the primary URL then any fallback URLs.
# Prints the HTTP status code of every attempt and, on success, the final size
# and elapsed time. Removes the temp file if the server returns an HTML error
# page (e.g. a 404 captured as a page) so the notebook never imports a bogus
# file. Returns 0 on success, 1 if every candidate failed.
# Usage: download "<primary_url>" "<out_file>" ["<fallback_url1>" ...]
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
# Build candidate ReMap storage URLs for a TF, trying letter-case variants
# (ReMap occasionally serves a file under a different letter-case).
# ------------------------------------------------------------------------------
remap_tf_candidates() {
  local tf="$1"
  local lc; lc=$(printf '%s' "$tf" | tr '[:upper:]' '[:lower:]')
  local tc; tc=$(printf '%s' "$tf" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')
  echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$tf/remap2022_${tf}_all_macs2_hg38_v1_0.bed.gz"
  [ "$lc" != "$tf" ] && echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$lc/remap2022_${lc}_all_macs2_hg38_v1_0.bed.gz"
  [ "$tc" != "$tf" ] && echo "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/$tc/remap2022_${tc}_all_macs2_hg38_v1_0.bed.gz"
}

echo "=============================================="
echo " Downloading external tracks -> $OUTDIR"
echo "=============================================="

# ------------------------------------------------------------------------------
# B1. RepeatMasker hg38 (UCSC)
# ------------------------------------------------------------------------------
echo "[B1] RepeatMasker hg38 (UCSC rmsk)"
download "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" "rmsk.hg38.bed.gz"

# ------------------------------------------------------------------------------
# B2. GTEx eQTL SNPs
# ------------------------------------------------------------------------------
# rsIDs from the report intersect STR loci. A representative set of eQTL SNPs
# cited in Relatorio_STR_Final_Integral.pdf (tissue in parentheses).
# Full GTEx eQTL data: https://gtexportal.org/home/downloads/adult-gtex/qtl
echo "[B2] GTEx eQTL SNP list (report-derived rsIDs)"
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

# ------------------------------------------------------------------------------
# B3. CTCF ChIP-seq signal (ENCODE, brain - dorsolateral prefrontal cortex)
#     ENCFF910VLV (fold change over control, GRCh38, released)
# ------------------------------------------------------------------------------
echo "[B3] CTCF ChIP-seq (ENCODE brain, ENCFF910VLV)"
download "https://www.encodeproject.org/files/ENCFF910VLV/@@download/ENCFF910VLV.bigWig" "CTCF_ENCFF910VLV.bw"

# ------------------------------------------------------------------------------
# B4. DNase-seq accessibility (ENCODE, brain)
#     ENCFF110ZAI (read-depth normalized signal, GRCh38, released)
# ------------------------------------------------------------------------------
echo "[B4] DNase-seq accessibility (ENCODE brain, ENCFF110ZAI)"
download "https://www.encodeproject.org/files/ENCFF110ZAI/@@download/ENCFF110ZAI.bigWig" "DNase_brain.bw"

# ------------------------------------------------------------------------------
# B5. H3K27ac histone mark (ENCODE, brain - dorsolateral prefrontal cortex)
#     ENCFF542HTC (fold change over control, GRCh38, released)
# ------------------------------------------------------------------------------
echo "[B5] H3K27ac (ENCODE brain, ENCFF542HTC)"
download "https://www.encodeproject.org/files/ENCFF542HTC/@@download/ENCFF542HTC.bigWig" "H3K27ac_brain.bw"

# ------------------------------------------------------------------------------
# B6. H3K27me3 repressive mark (ENCODE, brain) - replaces retired Roadmap link
#     (Roadmap Egg2 server no longer serves the DNA-methylation bigwigs.)
# ------------------------------------------------------------------------------
echo "[B6] H3K27me3 (ENCODE brain, ENCFF722UXM)"
download "https://www.encodeproject.org/files/ENCFF722UXM/@@download/ENCFF722UXM.bigWig" "H3K27me3_brain.bw"

# ------------------------------------------------------------------------------
# B7. TF ChIP-seq peaks (ReMap 2022) - CTCF, REST, KLF9, GATA2, ZNF384, MNT
#     Per-TF files: https://remap.univ-amu.fr/target_page/<TF>:9606
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# B8. ReMap Density - aggregated TF binding density (hg38)
#     Served by the ReMap 2022 UCSC track hub (the /storage/remap2022/.../density
#     path used to return 404).
# ------------------------------------------------------------------------------
echo "[B8] ReMap 2022 density (hg38)"
download "https://remap.univ-amu.fr/storage/public/hubReMap2022/hg38/bigBed/hg38.bw" "remap2022_density_hg38.bw"

# ------------------------------------------------------------------------------
# B9. ClinGen dosage sensitivity (haploinsufficiency / triplosensitivity)
#     Gene curation list (GRCh38) from ftp.clinicalgenome.org
# ------------------------------------------------------------------------------
echo "[B9] ClinGen dosage sensitivity (gene curation, GRCh38)"
download "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv" "ClinGen_gene_curation_list_GRCh38.tsv"

# ------------------------------------------------------------------------------
# B10. JARVIS (Junk Annotation genome-wide Residual Variation Intolerance Score)
#      Constraint score for non-coding regions (UCSC gbdb, hg38).
#      CAUTION: 21 GB bigWig. The notebook queries it by window over the remote
#      URL (no local download required), so this step is OPTIONAL.
#      Set DOWNLOAD_LARGE=1 to fetch it anyway, e.g.:
#        DOWNLOAD_LARGE=1 bash 7.4.3.1_download_external_tracks.sh
# ------------------------------------------------------------------------------
if [ "${DOWNLOAD_LARGE:-0}" = "1" ]; then
  echo "[B10] JARVIS non-coding constraint score (hg38, 21GB bigWig)"
  download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/jarvis/jarvis.bw" "jarvis.bw"
else
  echo "[B10] SKIP JARVIS bigWig (21 GB). Notebook queries remote URL by window."
fi

# ------------------------------------------------------------------------------
# B11. UK Biobank / deCODE Genetics Depletion Rank Score (hg38)
#      Constraint score: rank 0 (most depleted) to 100 (least depleted).
# ------------------------------------------------------------------------------
echo "[B11] UK Biobank / deCODE Depletion Rank Score (hg38)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/ukbDepletion/ukbDepletion.bw" "ukbDepletion.bw"

# ------------------------------------------------------------------------------
# B12. ReMap Atlas of Regulatory Regions - aggregated density (hg38, UCSC)
#      Distinct from the remap2022 density bigWig already downloaded in B8.
# ------------------------------------------------------------------------------
echo "[B12] ReMap Atlas density (hg38, UCSC reMapDensity2022)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMapDensity2022.bw" "remap_density_ucsc.bw"

# ------------------------------------------------------------------------------
# B13. JASPAR CORE TFBS predictions (hg38) - bigBed (billions of items)
#      Must be queried by window with bigBedToBed (installed via micromamba).
# ------------------------------------------------------------------------------
echo "[B13] JASPAR CORE TFBS (hg38, JASPAR2024.bb)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2024.bb" "JASPAR2024.bb"

# ------------------------------------------------------------------------------
# B14. DECIPHER common CNVs (hg38) - dbVar curated, AF>=0.01
#      bigBed; queried by window with bigBedToBed.
# ------------------------------------------------------------------------------
echo "[B14] DECIPHER common CNVs (hg38, common_decipher.bb)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/dbVar/common_decipher.bb" "decipher_common.bb"

# ------------------------------------------------------------------------------
# B15. TRExplorer tandem repeat catalog (hg38, UCSC strVar)
#      bigBed; queried by window with bigBedToBed.
# ------------------------------------------------------------------------------
echo "[B15] TRExplorer tandem repeat catalog (hg38, trexplorer.bb)"
download "https://hgdownload.soe.ucsc.edu/gbdb/hg38/strVar/trexplorer.bb" "trexplorer.bb"

# ------------------------------------------------------------------------------
# B16. Extra ReMap TF ChIP-seq peaks (beyond the original 6)
#      NACC2, FOXB1, FOXK1, PATZ1, RAD21, TCF4, MITF, USF2, RAD51, GATA1
# ------------------------------------------------------------------------------
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
# Best-effort: a missing TF peak (e.g. NACC2/FOXB1 absent from ReMap 2022 hg38)
# must NOT abort the whole run. Failures are collected and reported at the end.
declare -a FAILED_TFS=()
for tf in "${!EXTRA_TF_FILES[@]}"; do
  echo "   $tf"
  mapfile -t cands < <(remap_tf_candidates "$tf")
  if download "${cands[0]}" "${EXTRA_TF_FILES[$tf]}" "${cands[@]:1}"; then
    :
  else
    echo "   [WARN] $tf unavailable (HTTP 404 / not on ReMap 2022 hg38). Kept in list; will be handled separately."
    FAILED_TFS+=("$tf")
  fi
done

# ------------------------------------------------------------------------------
# B17. ENCODE wgEncodeReg4TfChip CTCF ChIP-seq (wgEncodeReg4TfChip_ENCFF341CQE)
#      BigWig from the ENCODE Portal accession ENCFF341CQE.
# ------------------------------------------------------------------------------
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
  echo "      investigate alternative sources or remove them from the list.)"
fi
ls -lh
echo "=============================================="
