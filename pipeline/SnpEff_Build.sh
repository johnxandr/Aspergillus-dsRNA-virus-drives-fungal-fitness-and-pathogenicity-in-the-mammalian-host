#!/bin/bash
#SBATCH --job-name=snpeff_reference
#SBATCH --output=snpeff_reference.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00

# Load Singularity module
module load singularity

# Set paths
SINGULARITY_IMAGE_URL="https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_0"
BASE_DIR="/sci/home/john.adeoye01/Desktop/mycoviruses/"
REF_GENOME_DIR="${BASE_DIR}ref_genome/"
GFF_FILE="${REF_GENOME_DIR}Afum293.gff"
SNPEFF_DATA_DIR="${BASE_DIR}snpEff_data/"  # Specify your desired directory under /sci/home

# Create snpEff_data directory if it doesn't exist
mkdir -p $SNPEFF_DATA_DIR

# Download Singularity image (if not already downloaded)
singularity pull $SINGULARITY_IMAGE_URL

# Run snpEff to build the reference genome
singularity run $SINGULARITY_IMAGE_URL snpEff build -gff3 -v Af293 -datadir $SNPEFF_DATA_DIR $GFF_FILE

