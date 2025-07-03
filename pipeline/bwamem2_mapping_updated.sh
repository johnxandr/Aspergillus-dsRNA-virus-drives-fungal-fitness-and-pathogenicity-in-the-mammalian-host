#!/bin/bash
#SBATCH --job-name=bwa_mem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=bwa_mem2_%A_%a.out
#SBATCH --error=bwa_mem2_%A_%a.err

module load singularity

# Define input and output directories
input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/FastP"
output_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/mapped"
individual_bam_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/individual_bam"

# Create output and individual BAM directories if they don't exist
mkdir -p $output_dir
mkdir -p $individual_bam_dir

# Load samtools module
module load samtools

# Define the reference genome path (replace with your actual reference genome path)
reference_genome="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"

# Index the reference genome with bwa
singularity run https://depot.galaxyproject.org/singularity/bwa-mem2%3A2.2.1--hd03093a_5 \
    bwa-mem2 index $reference_genome

# Loop through each sample and run bwa-mem2
for sample_dir in $input_dir/*; do
    sample_name=$(basename $sample_dir)
    
    # Create a folder for each sample in the mapped directory
    sample_output_dir="$output_dir/$sample_name"
    mkdir -p $sample_output_dir

    # Define paths for input reads
    read1="$sample_dir/${sample_name}_R1_trimmed.fastq.gz"
    read2="$sample_dir/${sample_name}_R2_trimmed.fastq.gz"

    # Define paths for output files
    output_sam="$sample_output_dir/${sample_name}_aligned.sam"
    output_bam="$sample_output_dir/${sample_name}_sorted.bam"
    output_bai="$sample_output_dir/${sample_name}_sorted.bam.bai"
    output_flagstat="$sample_output_dir/${sample_name}_flagstat.txt"
    individual_bam="$individual_bam_dir/${sample_name}_sorted.bam"

    # Run bwa-mem2
   singularity run https://depot.galaxyproject.org/singularity/bwa-mem2%3A2.2.1--hd03093a_5 \
    bwa-mem2 mem -t $SLURM_NTASKS -R "@RG\tID:$sample_name\tSM:$sample_name\tLB:library1\tPL:illumina" $reference_genome $read1 $read2 > $output_sam

    # Convert to BAM and sort using samtools
    samtools view -Sb $output_sam | samtools sort -o $output_bam -

    # Index the sorted BAM file
    samtools index $output_bam

    # Move the index file to the correct location
    mv ${output_bam}.bai $output_bai

    # Copy the individual BAM file to the separate directory
    cp $output_bam $individual_bam

    # Generate SAMtools flagstat
    samtools flagstat $output_bam > $output_flagstat
done


