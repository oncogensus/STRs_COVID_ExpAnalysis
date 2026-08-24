#!/bin/bash
# ==============================================================================
# 7.4.3.1_debug_nacc2_foxb1.sh
# ==============================================================================
# Debug script to search for alternative ChIP-seq sources for NACC2 and FOXB1
# in hg38/GRCh38, since they are absent from ReMap 2022.
#
# Sources checked:
#   - ENCODE Portal (via API)
#   - Cistrome DB
#   - ReMap target pages (direct)
#   - UCSC public track hubs
#   - GEO/SRA (if needed)
#
# Usage:
#   bash 7.4.3.1_debug_nacc2_foxb1.sh
#
# Output: debug log + suggested URLs for main script
# ==============================================================================
set -euo pipefail

OUTDIR="external_tracks"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOGFILE="debug_nacc2_foxb1_$(date '+%Y%m%d_%H%M%S').log"
exec > >(
    while IFS= read -r line; do
        echo "$(date '+%Y-%m-%d %H:%M:%S') $line"
    done | tee -a "$LOGFILE"
) 2>&1

echo "=============================================="
echo " Debug: Alternative sources for NACC2/FOXB1"
echo " Log: $LOGFILE"
echo "=============================================="

# TFs to search for
declare -A TFS=(
    [NACC2]="NACC2"
    [FOXB1]="FOXB1"
)

# Base URLs to try
declare -A BASE_URLS=(
    [remap_storage]="https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF"
    [encode]="https://www.encodeproject.org"
    [cistrome]="http://cistrome.org/db"
    [ucsc_hub]="https://hgdownload.soe.ucsc.edu/hubs"
)

# Search function for a single TF
search_tf() {
    local tf="$1"
    local lc_tf; lc_tf=$(printf '%s' "$tf" | tr '[:upper:]' '[:lower:]')
    local tc_tf; tc_tf=$(printf '%s' "$tf" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')

    echo ""
    echo "=== Searching for $tf ==="

    # 1. Try ReMap with different casings (already done in main script, but re-verify)
    echo "1. ReMap storage (casings):"
    for variant in "$tf" "$lc_tf" "$tc_tf"; do
        url="${BASE_URLS[remap_storage]}/$variant/remap2022_${variant}_all_macs2_hg38_v1_0.bed.gz"
        http_code=$(curl -sS -o /dev/null -w '%{http_code}' -I "$url" --max-time 10 2>/dev/null || echo "000")
        echo "   $variant: HTTP $http_code"
        if [ "$http_code" = "200" ]; then
            echo "   >>> FOUND: $url"
            return 0
        fi
    done

    # 2. Check ReMap target page for this TF
    echo "2. ReMap target page:"
    target_url="https://remap.univ-amu.fr/target_page/${tf}:9606"
    http_code=$(curl -sS -o /dev/null -w '%{http_code}' -I "$target_url" --max-time 10 2>/dev/null || echo "000")
    echo "   $target_url: HTTP $http_code"
    if [ "$http_code" = "200" ]; then
        echo "   >>> Target page exists, check for download links"
        # Could parse HTML for actual file links
    fi

    # 3. ENCODE search via API
    echo "3. ENCODE Portal search:"
    search_url="https://www.encodeproject.org/search/?type=Experiment&target.title=${tf}&assembly=GRCh38&status=released&format=json"
    response=$(curl -sS "$search_url" --max-time 15 2>/dev/null || echo '{"@graph":[]}')
    count=$(echo "$response" | jq '.["@graph"] | length' 2>/dev/null || echo 0)
    echo "   ENCODE experiments found: $count"
    if [ "$count" -gt 0 ]; then
        echo "$response" | jq -r '.["@graph"][] | "\(.accession) \(.target.title) \(.biosample_ontology.term_name // "N/A")"' 2>/dev/null | head -10 | while IFS= read -r line; do
            echo "   $line"
        done
        # Could extract bigWig/bed.gz file URLs from experiments
    fi

    # 4. Cistrome DB check (simple HTTP check)
    echo "4. Cistrome DB (simple check):"
    cistrome_url="http://cistrome.org/db/#${tf}"
    http_code=$(curl -sS -o /dev/null -w '%{http_code}' -I "$cistrome_url" --max-time 10 2>/dev/null || echo "000")
    echo "   $cistrome_url: HTTP $http_code"

    # 5. UCSC public track hubs check
    echo "5. UCSC public hubs (general check):"
    for hub in "encode" "remap" "roadmap"; do
        hub_url="https://hgdownload.soe.ucsc.edu/hubs/${hub}/hub.txt"
        http_code=$(curl -sS -o /dev/null -w '%{http_code}' -I "$hub_url" --max-time 10 2>/dev/null || echo "000")
        echo "   $hub hub: HTTP $http_code"
    done

    echo ""
    return 1
}

# Main search loop
for tf in "${!TFS[@]}"; do
    search_tf "$tf"
done

echo "=============================================="
echo " Search complete. Log: $LOGFILE"
echo "=============================================="

# Summary of findings
echo ""
echo "=== SUMMARY ==="
echo "If any source returned HTTP 200, note the URL above and add to main script."
echo "If all returned 404, the TFs are genuinely absent from those databases."
echo ""
echo "Next steps if alternatives found:"
echo "  1. Add URLs to 7.4.3.1_download_external_tracks.sh EXTRA_TF_FILES"
echo "  2. Or create separate download step for these TFs"
echo ""
echo "Log saved to: $LOGFILE"