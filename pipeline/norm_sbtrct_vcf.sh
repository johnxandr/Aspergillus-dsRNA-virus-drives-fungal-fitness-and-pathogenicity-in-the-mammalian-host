#!/bin/bash
#SBATCH --job-name=bcftools_subtract
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=bcftools_subtract.out
#SBATCH --error=bcftools_subtract.err

# Load bcftools module
module load bcftools

# Set input and output directories
input_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls"
output_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/subtracted"

# Create output directory if not exists
mkdir -p $output_dir

# Normalize and subtract for each sample
for sample_dir in $input_dir/*; do
  if [ -d "$sample_dir" ]; then
    sample_name=$(basename "$sample_dir")
    input_vcf="$sample_dir/${sample_name}_raw_variants.g.vcf.gz"

    # Perform normalization and subtraction using bcftools
    bcftools norm -Ou -m -any $input_vcf | \
    bcftools isec -c none -Oz -p "$output_dir/$sample_name" \
      - $input_dir/anc_gmm/anc_gmm_raw_variants.g.vcf.gz

    # Index the resulting VCF file
    # bcftools index "$output_dir/$sample_name/$sample_name.0000.vcf.gz"
  fi
done
