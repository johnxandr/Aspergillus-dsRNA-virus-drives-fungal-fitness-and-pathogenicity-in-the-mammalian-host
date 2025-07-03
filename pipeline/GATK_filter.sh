#!/bin/bash
#SBATCH --job-name=gatk_filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=gatk_filter.out
#SBATCH --error=gatk_filter.err

module load singularity
module load bcftools

# Set the input and output directories
reference_genome="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls"
output_root="/sci/home/john.adeoye01/Desktop/mycoviruses/filtered_gatk"

# Load the GATK Singularity image
gatk_image="https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"

# Iterate over sample directories
for sample_dir in "$input_dir"/*; do
    sample_name=$(basename "$sample_dir")
    output_dir="$output_root/$sample_name"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Task 1: Select the SNPs with GATK command
    snp_output="$output_dir/${sample_name}_snp_output.vcf"
    singularity exec "$gatk_image" gatk SelectVariants \
        -R "$reference_genome" \
        -V "$sample_dir/${sample_name}_raw_variants.g.vcf" \
        --select-type SNP \
        -O "$snp_output"

    # Task 2: Filter the SNPs
    snp_filter="--filter-expression 'QD < 3.0' --filter-expression 'MQ < 40.0' --filter-expression 'FS > 60.0' --filter-expression 'MQRankSum < -10.0' --filter-expression 'ReadPosRankSum < -5.0' --filter-name lowQC"
    filtered_snp_output="$output_dir/${sample_name}_filtered_snps.vcf"
    singularity exec "$gatk_image" gatk VariantFiltration \
        -V "$snp_output" \
        -O "$filtered_snp_output" \
        --filter "$snp_filter"

    # Task 3: Select the INDELs with the GATK command
    indel_output="$output_dir/${sample_name}_indel_output.vcf"
    singularity exec "$gatk_image" gatk SelectVariants \
        -R "$reference_genome" \
        -V "$sample_dir/${sample_name}_raw_variants.g.vcf" \
        --select-type INDEL \
        -O "$indel_output"

    # Task 4: Filter the INDELs
    indel_filter="--filter-expression 'QD < 3.0' --filter-expression 'FS > 200.0' --filter-expression 'ReadPosRankSum < -20.0' --filter-name lowQC"
    filtered_indel_output="$output_dir/${sample_name}_filtered_indels.vcf"
    singularity exec "$gatk_image" gatk VariantFiltration \
        -V "$indel_output" \
        -O "$filtered_indel_output" \
        --filter "$indel_filter"

    # Task 5: Merge SNP and INDEL VCFs using bcftools
    merged_output="$output_dir/${sample_name}_merged_snps_indels.vcf"
    bcftools merge "$filtered_snp_output" "$filtered_indel_output" -O v -o "$merged_output"
done

