#!/bin/bash
#SBATCH --job-name=platypus_variant_calling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8  # Adjust according to your system's specifications
#SBATCH --mem=16G  # Adjust according to your system's specifications
#SBATCH --time=24:00:00  # Adjust according to your needs
#SBATCH --output=platypus_variant_calling_%j.out
#SBATCH --error=platypus_variant_calling_%j.err

# Load required module
module load singularity

# Set paths
REFERENCE_GENOME="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
INPUT_DIR="/sci/home/john.adeoye01/Desktop/mycoviruses/Platypus"
OUTPUT_DIR="/sci/home/john.adeoye01/Desktop/mycoviruses/platypus_variants_2"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Singularity exec with URL
PLATYPUS_IMAGE_URL="https://depot.galaxyproject.org/singularity/platypus-variant%3A0.8.1.2--py27hb698ca4_5"

# Quality filter options (Adjust these values based on your criteria)
QUAL_FILTER_THRESHOLD=30
DP_FILTER_THRESHOLD=10
AF_FILTER_THRESHOLD=0.1

# Iterate through each sample and perform variant calling
for BAM_FILE in "$INPUT_DIR"/*/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" | cut -d '_' -f 1)
    OUTPUT_VCF="$OUTPUT_DIR/${SAMPLE_NAME}_platypus.vcf"
    
    singularity exec --bind "$REFERENCE_GENOME":"$REFERENCE_GENOME" "$PLATYPUS_IMAGE_URL" platypus callVariants \
        --bamFiles="$BAM_FILE" \
        --refFile="$REFERENCE_GENOME" \
        --output="$OUTPUT_VCF" \
        --nCPU=8  # Adjust according to your system's specifications

    # Add Platypus filtering options
    singularity exec "$PLATYPUS_IMAGE_URL" platypus filterVariants \
        --vcf "$OUTPUT_VCF" \
        --minQuality="$QUAL_FILTER_THRESHOLD" \
        --minDepth="$DP_FILTER_THRESHOLD" \
        --minVarFreq="$AF_FILTER_THRESHOLD" \
        --output "$OUTPUT_VCF.filtered"

    echo "Variant calling and filtering completed for $SAMPLE_NAME. Output VCF: $OUTPUT_VCF.filtered"
done

echo "Variant calling and filtering completed for all samples!"
