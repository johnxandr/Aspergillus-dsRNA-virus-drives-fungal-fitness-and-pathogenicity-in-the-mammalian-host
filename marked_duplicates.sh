#!/bin/bash
#SBATCH --job-name=mark_duplicates
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=mark_duplicates_%A.out
#SBATCH --error=mark_duplicates_%A.err

module load singularity
module load samtools

# Define input and output directories
input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/mapped"
output_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/marked_duplicates"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Load GATK4 Singularity image
gatk_singularity_image="https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"

# Loop through each sample and mark duplicates
for sample_dir in $input_dir/*; do
    sample_name=$(basename $sample_dir)

    # Define paths for input and output files
    input_bam="$sample_dir/${sample_name}_sorted.bam"
    output_bam="$output_dir/${sample_name}_dedup.bam"
    metrics_file="$output_dir/${sample_name}_dedup_metrics.txt"

    # Add read groups using Picard
    # singularity exec $gatk_singularity_image gatk MarkDuplicatesSpark \
       #  I=$input_bam \
       #  O=$output_bam \
       #  RGID=1 \
       #  RGLB=library \
       #  RGPL=platform \
       #  RGSM=$sample_name             

   # Use GATK4's MarkDuplicatesSpark within Singularity
    singularity exec $gatk_singularity_image gatk MarkDuplicatesSpark \
       --input $input_bam \
       --output $output_bam \
       --metrics-file $metrics_file \
       --remove-sequencing-duplicates false  # Set to true if you want to remove duplicates, false to just mark them \
       # --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g"  # Adjust memory as needed
   # Index the marked duplicates BAM file using samtools within Singularity
    singularity exec $gatk_singularity_image module load samtools && samtools index $output_bam
done

