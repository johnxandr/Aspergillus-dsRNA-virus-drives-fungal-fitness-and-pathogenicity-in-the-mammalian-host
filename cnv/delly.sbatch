#!/bin/bash
# -----------------------------------------------------------------------------
# Description:        DELLY analysis
# Author:             Maria Laura Fabre (https://github.com/laurafabre)
# Date Created:       <2024-03-08>
# Last Modified:      <2024-03-08>
# Usage:              bash delly.sbatch 
# Notes:             It is needed to hcreate first the mappability map
# -----------------------------------------------------------------------------
# Get the directory path of the script
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUTPUT_DIR="$SCRIPT_DIR/results/variant_calling/delly/"
REF_DIR="$SCRIPT_DIR/results/variant_calling/delly/ref_genome"
MAP_DIR="$SCRIPT_DIR/results/variant_calling/delly/mappability"
BAM_DIR="freec"


# Creation Mappability map
# copy reference genome to delly folder
mkdir -p "$REF_DIR"
mkdir -p "$MAP_DIR"
cp "$SCRIPT_DIR/ref_genome/Afum293_genome.fasta" "$REF_DIR/"

# Dicey Chop
cd "$OUTPUT_DIR"
dicey chop "$REF_DIR/Afum293_genome.fasta"

# Move the auto generated reads
mv read1.fq.gz read2.fq.gz "$MAP_DIR/"
chmod +x "$MAP_DIR/read*"

# Index reference
bwa index "$REF_DIR/Afum293_genome.fasta"

# Align reads and sort

bwa mem "$REF_DIR/Afum293_genome.fasta" "$MAP_DIR/read1.fq.gz" "$MAP_DIR/read2.fq.gz" | samtools sort -@ 8 -o "$REF_DIR/Afum293_genome.fasta.srt.bam" 

# Index sorted BAM
samtools index "$REF_DIR/Afum293_genome.fasta.srt.bam"

# Calculate mappability
dicey mappability2 "$REF_DIR/Afum293_genome.fasta.srt.bam" 
mv "$OUTPUT_DIR/map.fa.gz" "$MAP_DIR/"
gunzip -c "$MAP_DIR/map.fa.gz" > "$MAP_DIR/map.fa"
bgzip "$MAP_DIR/map.fa"
samtools faidx "$MAP_DIR/map.fa.gz"

# Delly CNV analysis

# Find all .bam files in subfolders and perform Delly CNV analysis
find "$BAM_DIR" -type f -name "*.bam" | while read -r BAM_FILE; do
    # Extract sample name from the path
    SAMPLE_NAME=$(basename "${BAM_FILE%.*}") # Remove the .bam extension
    
    # Perform CNV analysis with Delly for each .bam file found
    delly cnv -o  "$OUTPUT_DIR/$SAMPLE_NAME.bcf" --ploidy 1 -g "$REF_DIR/Afum293_genome.fasta" -m "$MAP_DIR/map.fa.gz" "$BAM_FILE"
    
    # Convert BCF to VCF
    bcftools convert -O v -o "$OUTPUT_DIR/$SAMPLE_NAME.vcf" "$OUTPUT_DIR/$SAMPLE_NAME.bcf"
done
