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
#   B3. CTCF ChIP-seq (ENCODE) (wgEncodeReg4TfChip_ENCFF341CQE)
#   B4. DNase-seq              (ENCODE)
#   B5. H3K27ac                (ENCODE / Roadmap)
#   B6. DNA methylation        (ENCODE/Roadmap, brain)
#   B7. TF ChIP-seq (ReMap)    (NACC2, FOXB1, KLF9, GATA2, ZNF384, MNT)
#   B8. ReMap Density          (ReMap 2022)
#   B9. DECIPHER               (dosage sensitivity)
#
# Usage:
#   bash 7.4.3.1_download_external_tracks.sh
#
# Output: all files saved under ./external_tracks/
# ==============================================================================
set -euo pipefail

OUTDIR="external_tracks"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=============================================="
echo " Downloading external tracks -> $OUTDIR"
echo "=============================================="

# ------------------------------------------------------------------------------
# B1. RepeatMasker hg38 (UCSC)
# ------------------------------------------------------------------------------
echo "[B1] RepeatMasker hg38 (UCSC rmsk)"
if [ ! -f rmsk.hg38.bed.gz ]; then
  curl -sSL -o rmsk.hg38.bed.gz \
    "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
fi

# ------------------------------------------------------------------------------
# B2. GTEx eQTL SNPs
# ------------------------------------------------------------------------------
# rsIDs from the report intersect STR loci. A representative set of eQTL SNPs
# cited in Relatorio_STR_Final_Integral.pdf (tissue in parentheses).
# Full GTEx eQTL data: https://gtexportal.org/home/downloads/adult-gtex/qtl
echo "[B2] GTEx eQTL SNP list (report-derived rsIDs)"
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
echo "   -> saved gtex_eqTL_rsids.tsv ($(wc -l < gtex_eqTL_rsids.tsv) lines)"

# ------------------------------------------------------------------------------
# B3. CTCF ChIP-seq signal (ENCODE, wgEncodeReg4TfChip_ENCFF341CQE)
# ------------------------------------------------------------------------------
echo "[B3] CTCF ChIP-seq (ENCODE wgEncodeReg4TfChip_ENCFF341CQE)"
if [ ! -f CTCF_ENCFF341CQE.bw ]; then
  curl -sSL -o CTCF_ENCFF341CQE.bw \
    "https://www.encodeproject.org/files/ENCFF341CQE/@@download/ENCFF341CQE.bigWig"
fi

# ------------------------------------------------------------------------------
# B4. DNase-seq accessibility (ENCODE, brain) - representative track
# ------------------------------------------------------------------------------
echo "[B4] DNase-seq accessibility (ENCODE, brain)"
if [ ! -f DNase_brain.bw ]; then
  curl -sSL -o DNase_brain.bw \
    "https://www.encodeproject.org/files/ENCFF284JLI/@@download/ENCFF284JLI.bigWig"
fi

# ------------------------------------------------------------------------------
# B5. H3K27ac histone mark (ENCODE/Roadmap, brain) - representative track
# ------------------------------------------------------------------------------
echo "[B5] H3K27ac (ENCODE, brain)"
if [ ! -f H3K27ac_brain.bw ]; then
  curl -sSL -o H3K27ac_brain.bw \
    "https://www.encodeproject.org/files/ENCFF744IFL/@@download/ENCFF744IFL.bigWig"
fi

# ------------------------------------------------------------------------------
# B6. DNA methylation (Roadmap, brain) - representative track
# ------------------------------------------------------------------------------
echo "[B6] DNA methylation (Roadmap brain)"
if [ ! -f Methylation_brain.bed.gz ]; then
  curl -sSL -o Methylation_brain.bed.gz \
    "https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/MethylHMR/FractionalMethylation_bigwig/E002-H3K4me1.bigwig"
fi

# ------------------------------------------------------------------------------
# B7. TF ChIP-seq peaks (ReMap 2022) - NACC2, FOXB1, KLF9, GATA2, ZNF384, MNT
#     ReMap peak files: https://remap.univ-amu.fr/
# ------------------------------------------------------------------------------
echo "[B7] ReMap TF ChIP-seq peaks"
declare -A TF_FILES=(
  [NACC2]="remap2022_nacc2_all_macs2_hg38_v1_0.bed.gz"
  [FOXB1]="remap2022_foxb1_all_macs2_hg38_v1_0.bed.gz"
  [KLF9]="remap2022_klf9_all_macs2_hg38_v1_0.bed.gz"
  [GATA2]="remap2022_gata2_all_macs2_hg38_v1_0.bed.gz"
  [ZNF384]="remap2022_znf384_all_macs2_hg38_v1_0.bed.gz"
  [MNT]="remap2022_mnt_all_macs2_hg38_v1_0.bed.gz"
)
for tf in "${!TF_FILES[@]}"; do
  if [ ! -f "${TF_FILES[$tf]}" ]; then
    echo "   downloading $tf"
    curl -sSL -o "${TF_FILES[$tf]}" \
      "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/ALL/${TF_FILES[$tf]}"
  fi
done

# ------------------------------------------------------------------------------
# B8. ReMap Density - aggregated TF binding density (hg38)
# ------------------------------------------------------------------------------
echo "[B8] ReMap 2022 density (hg38)"
if [ ! -f remap2022_density_hg38.bw ]; then
  curl -sSL -o remap2022_density_hg38.bw \
    "https://remap.univ-amu.fr/storage/remap2022/hg38/density/remap2022_density_hg38_v1_0.bw"
fi

# ------------------------------------------------------------------------------
# B9. DECIPHER dosage sensitivity (haploinsufficiency / triplosensitivity)
#     Public ClinGen dosage sensitivity download
# ------------------------------------------------------------------------------
echo "[B9] DECIPHER/ClinGen dosage sensitivity"
if [ ! -f ClinGen_dosage_sensitivity.csv ]; then
  curl -sSL -o ClinGen_dosage_sensitivity.csv \
    "https://ftp.clinicalgenome.org/ClinGen_dosage_sensitivity.csv"
fi

echo ""
echo "=============================================="
echo " DONE. Tracks available in: $(pwd)"
ls -lh
echo "=============================================="
