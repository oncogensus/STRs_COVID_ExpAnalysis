#!/bin/bash

# Configuration variables
PREFIX="*"
PREFIX_PATH="/storage/users/tulio/Projeto_Luy_COVID/results/recal"
STRLING="/storage2/matheusbomfim/projects/micromamba/envs/str/bin/strling"
REF_FA="hg38.fa"
GENOME_STR="hg38.fa.str"
JOINT_DIR="cbgm_output"

# Create a directory if necessary
mkdir -p $JOINT_DIR

echo "Starting file extraction ${PREFIX}*.bam in the directory: $PREFIX_PATH..."

# Enable nullglob so that the array is empty if there are no files
shopt -s nullglob
files=(${PREFIX_PATH}/${PREFIX}*.bam)

if [ ${#files[@]} -eq 0 ]; then
    echo "No files ${PREFIX}*.bam found in $PREFIX_PATH. Aborting."
    exit 1
fi

for i in "${files[@]}"; do
    basefile=$(basename "$i")
    echo "Extracting STRs from the file $basefile..."
    $STRLING extract -f $REF_FA -g $GENOME_STR -v "$i" $JOINT_DIR/"$basefile".str.bin
    if [ $? -ne 0 ]; then
        echo "Error extracting the file $basefile. Aborting."
        exit 1
    fi
done

echo "Extraction successfully completed."

echo "Starting merge of .str.bin files..."

$STRLING merge -f $REF_FA -v -o $JOINT_DIR/strling $JOINT_DIR/*.str.bin
if [ $? -ne 0 ]; then
    echo "Merge error. Aborting."
    exit 1
fi

echo "Merge successfully completed."

echo "Starting STR calling on files ${PREFIX}*.bam..."

for i in "${files[@]}"; do
    basefile=$(basename "$i")
    echo "Calling STRs on the file $basefile..."
    $STRLING call -v -b $JOINT_DIR/strling-bounds.txt -o $JOINT_DIR/"$basefile" "$i" $JOINT_DIR/"$basefile".str.bin
    if [ $? -ne 0 ]; then
        echo "Error calling STRs on the file $basefile. Aborting."
        exit 1
    fi
done

echo "Process completed successfully."
