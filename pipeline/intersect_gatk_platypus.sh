#!/bin/bash
#SBATCH --job-name=mycoviruses_processing
#SBATCH --output=mycoviruses_processing_%j.out
#SBATCH --error=mycoviruses_processing_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load bcftools

# Define root directory and reference genome
root_directory='/sci/home/john.adeoye01/Desktop/mycoviruses/intersec_platypus_gatk'
reference_genome='/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta'

# Function to compress and index VCF files
function compress_and_index_vcf {
    if [ -f "${1}.gz" ]; then
        rm "${1}.gz"  # Remove existing file before compression
    fi
    bcftools view -Oz -o "${1}.gz" "$1"
    bcftools index "${1}.gz"
}

# 2. Compress and index GATK files
find "$root_directory/isec-gatk" -type f -name '*.sbkgd.vcf' ! -name '*_sbkgd.vcf' | while read -r gatk_file; do
    compress_and_index_vcf "$gatk_file"
done

# 3. Compress and index Platypus files
find "$root_directory/isec_platypus" -type f -name '*_sbkgd_platypus.vcf' | while read -r platypus_file; do
    compress_and_index_vcf "$platypus_file"
done

# 4. Perform intersection for corresponding samples of GATK and Platypus files
intersection_dir="$root_directory/intersection_output"

# Iterate over GATK files
find "$root_directory/isec-gatk" -type f -name '*.sbkgd.vcf' ! -name '*_sbkgd.vcf' | while read -r gatk_file; do
    subdir=$(dirname "$gatk_file" | sed "s|$root_directory/isec-gatk/||")
    platypus_file="$root_directory/isec_platypus/$subdir/$(basename "$gatk_file" | sed 's/.sbkgd.vcf/_sbkgd_platypus.vcf.gz/')"
    
    # Perform intersection for the pair
    bcftools isec -p "$intersection_dir/$subdir" -n=2 -w1 "$gatk_file.gz" "$platypus_file"
done

# 5. Rename specific output files for each sample
find "$intersection_dir" -type d -mindepth 1 -maxdepth 1 | while read -r sample_dir; do
    sample_name=$(basename "$sample_dir")
    mv "$sample_dir/0000.vcf" "$root_directory/${sample_name}_intersect.vcf"
    mv "$sample_dir/sites.txt" "$root_directory/${sample_name}_sites.txt"
done

