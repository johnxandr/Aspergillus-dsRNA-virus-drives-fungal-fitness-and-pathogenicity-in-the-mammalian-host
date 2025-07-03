#!/bin/bash
#SBATCH --job-name=gatk_bqsr
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00  # Adjust based on your data and resources
#SBATCH --output=gatk_bqsr.out
#SBATCH --error=gatk_bqsr.err

# Load required modules
module load singularity

# Set the paths and variables
gatk_singularity_image="https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
output_recalibrated_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/recalibrated/"
known_sites="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Af293_SNPs_fungidb_summary.sorted.removedNoise.vcf.gz"
reference_genome="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
reference_dict="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.dict"

# Create the Recalibrated directory if it doesn't exist
mkdir -p $output_recalibrated_dir

# Download the Singularity container
singularity pull $gatk_singularity_image

# Loop through only _dedup.bam files
for input_bam in /sci/home/john.adeoye01/Desktop/mycoviruses/marked_duplicates/*_dedup.bam; do
    # Extract sample name from the BAM file
    sample_name=$(basename "$input_bam" "_dedup.bam")

    # Create a directory for each sample if it doesn't exist
    mkdir -p $output_recalibrated_dir/$sample_name

    # Use a temporary file for known-sites content
    tmp_known_sites="/tmp/known_sites_tmp.vcf"
    zcat $known_sites > $tmp_known_sites

    # Index the VCF file
    singularity exec $gatk_singularity_image gatk IndexFeatureFile -I $tmp_known_sites

    # Build the model for GATK BaseRecalibrator with known SNPs
    singularity exec $gatk_singularity_image gatk --java-options "-Xmx4g" BaseRecalibrator \
        -I $input_bam \
        -R $reference_genome \
        --known-sites $tmp_known_sites \
        -O $output_recalibrated_dir/$sample_name/${sample_name}_recal_data.table \
        --sequence-dictionary $reference_dict

    # Run GATK ApplyBQSR
    singularity exec $gatk_singularity_image gatk --java-options "-Xmx4g" ApplyBQSR \
        -I $input_bam \
        -R $reference_genome \
        --bqsr-recal-file $output_recalibrated_dir/$sample_name/${sample_name}_recal_data.table \
        -O $output_recalibrated_dir/$sample_name/${sample_name}_recal.bam \
        --sequence-dictionary $reference_dict

    # Remove the temporary file
    rm $tmp_known_sites
done
