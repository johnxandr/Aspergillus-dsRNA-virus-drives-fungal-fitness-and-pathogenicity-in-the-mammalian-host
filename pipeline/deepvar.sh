#!/bin/bash
#SBATCH --job-name=deepvariant_calling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4  # Adjust based on your system's resources
#SBATCH --mem=16G  # Adjust based on your system's resources
#SBATCH --time=12:00:00  # Adjust based on your expected run time
#SBATCH --output=deepvariant_calling_%A.out
#SBATCH --error=deepvariant_calling_%A.err

# Load required module
module load singularity

# Set your directories
reference_genome="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/Platypus"
output_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/deep_variants"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Download DeepVariant Singularity image with increased timeout (7200 seconds)
singularity pull --timeout 7200 --name deepvariant.sif docker://google/deepvariant:1.4.0--py36hf3e76ba_0

# Set the local path to the Singularity image
singularity_image="https://depot.galaxyproject.org/singularity/deepvariant%3A1.4.0--py36hf3e76ba_0"

# Loop through each BAM file in the input directory
for bam_file in $input_dir/*/*.bam; do
    # Extract sample name from the BAM file
    sample_name=$(basename "$bam_file" | cut -d '_' -f 1)

    # Define output VCF file name
    output_vcf="$output_dir/${sample_name}_deepvariant.vcf"

    # Run DeepVariant using Singularity
    singularity run $singularity_image \
        /opt/deepvariant-1.4.0/bin/run_deepvariant \
        --model_type=WGS \
        --ref=$reference_genome \
        --reads=$bam_file \
        --output_vcf=$output_vcf \
        --num_shards=4  # Adjust based on your system's resources
done

echo "DeepVariant variant calling completed."

