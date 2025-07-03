#!/bin/bash
#SBATCH --job-name=variant_filtering
#SBATCH --output=variant_filtering_%j.out
#SBATCH --error=variant_filtering_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

# Load Singularity module
module load singularity

# Set the path to your Singularity image
SINGULARITY_IMAGE_URL=https://depot.galaxyproject.org/singularity/platypus-variant%3A0.8.1.2--py27hb698ca4_5

# Set the input directory for VCF files
INPUT_DIR=/sci/home/john.adeoye01/Desktop/mycoviruses/platypus_variants

# Set the output directory for filtered VCF files
OUTPUT_DIR=/sci/home/john.adeoye01/Desktop/mycoviruses/platypus_variants/filtered

# Platypus filtering options
PLATYPUS_FILTER_OPTIONS="--minReads 5
                         --minVarFreq 0.1
                         --maxHeterozygosity 0.8
                         --maxIndelSize 10
                         --minIndelSize 1"

# Create output directory if not exists
mkdir -p $OUTPUT_DIR

# Iterate over VCF files in the input directory and apply Platypus filtering
for VCF_FILE in $INPUT_DIR/*.vcf; do
    BASENAME=$(basename $VCF_FILE .vcf)
    OUTPUT_VCF=$OUTPUT_DIR/${BASENAME}_filtered.vcf

    # Run Platypus within Singularity container with filtering options
    singularity exec $SINGULARITY_IMAGE_URL env LC_ALL=C platypus filterVariants $PLATYPUS_FILTER_OPTIONS \
        --vcf $VCF_FILE \
        --output $OUTPUT_VCF

    echo "Filtered VCF for $BASENAME created: $OUTPUT_VCF"
done

echo "Variant filtering job completed."

