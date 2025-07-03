#!/bin/bash
#SBATCH --job-name=subtract_background
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=subtract_background.out
#SBATCH --error=subtract_background.err

# Load required module
module load bcftools

# Define input directory for variants called by GATK
haplotypecaller_input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/filtered_gatk/filtered_gatk_snps"
background_file="/sci/home/john.adeoye01/Desktop/mycoviruses/filtered_gatk/filtered_gatk_snps/anc_gmm/anc_gmm_snp_output.vcf"

# Task 1: Normalize, compress, and index the background strain
background_norm="${background_file}.norm"
bcftools norm -m -any "${background_file}" -o "${background_norm}"
bcftools view -Oz -o "${background_norm}.gz" "${background_norm}"
bcftools index "${background_norm}.gz"

# Task 2: Normalize, compress with bgzip, and index the VCF files of the sample strains
sample_files=("RI" "VC-Cyclohexamide" "VC-Ribavirin")

for sample in "${sample_files[@]}"; do
    vcf_file="${haplotypecaller_input_dir}/${sample}/${sample}_snp_output.vcf"
    normalized_file="${vcf_file}.norm"
    
    bcftools norm -m -any "${vcf_file}" -o "${normalized_file}"
    bcftools view -Oz -o "${normalized_file}.gz" "${normalized_file}"
    bcftools index "${normalized_file}.gz"
done

# Task 3: Create isec folder to store the subtracted VCF files
for sample in "${sample_files[@]}"; do
    isec_output_dir="${haplotypecaller_input_dir}/isec/${sample}"
    mkdir -p "${isec_output_dir}"
done

# Task 4: Perform intersection with sample and background files
for sample in "${sample_files[@]}"; do
    vcf_file="${haplotypecaller_input_dir}/${sample}/${sample}_snp_output.vcf"
    normalized_file="${vcf_file}.norm.gz"
    
    bcftools isec --complement -p "${haplotypecaller_input_dir}/isec/${sample}" "${normalized_file}" "${background_norm}.gz"
done

# Task 5: Rename and move the output files
for sample in "${sample_files[@]}"; do
    isec_output_dir="${haplotypecaller_input_dir}/isec/${sample}"
    mv "${isec_output_dir}/0000.vcf" "${haplotypecaller_input_dir}/${sample}/${sample}.sbkgd.vcf"
    mv "${isec_output_dir}/sites.txt" "${haplotypecaller_input_dir}/${sample}/${sample}.sbkgd.txt"
done

