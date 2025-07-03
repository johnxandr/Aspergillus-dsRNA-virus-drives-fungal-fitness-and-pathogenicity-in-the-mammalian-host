#!/bin/bash
#SBATCH --job-name=fastp_preprocessing
#SBATCH --output=fastp_preprocessing.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00

# Load required module
module load singularity

# Set input and output base directories
input_base="/sci/home/john.adeoye01/Desktop/mycoviruses/rawReads"
output_base="/sci/home/john.adeoye01/Desktop/mycoviruses/FastP"

# Get sample names from subdirectories in the input base directory
samples=($(find "$input_base" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))

# Create output directories and run fastp for each sample
for sample in "${samples[@]}"; do
    input_dir="$input_base/$sample"
    output_dir="$output_base/$sample"

    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Unzip the FASTQ files
    #gzip "$input_dir/${sample}_R1_001.fastq"
    #gzip "$input_dir/${sample}_R2_001.fastq"

    # Run fastp using singularity container
    singularity run https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_2 \
      fastp -i "$input_dir/${sample}_R1_001.fastq.gz" -I "$input_dir/${sample}_R2_001.fastq.gz" \
      --out1 "$output_dir/${sample}_R1_trimmed.fastq.gz" --out2 "$output_dir/${sample}_R2_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 30 \
      --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      --overrepresentation_analysis --overrepresentation_sampling 20 \
      --low_complexity_filter --complexity_threshold 30 \
      --length_required 100 \
      --trim_poly_g --correction -f 5 -t 5 -F 5 -T 5 --trim_poly_x

    echo "FastP preprocessing for $sample complete."
done

echo "FastP preprocessing for all samples complete."

