#!/bin/bash
#SBATCH --job-name=gatk_variant_calling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00  # Adjust based on your data and resources
#SBATCH --output=gatk_variant_calling.out
#SBATCH --error=gatk_variant_calling.err

# Load required modules
module load singularity

# Set the paths and variables
gatk_singularity_image="https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
output_variant_dir="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/"
reference_genome="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
reference_dict="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.dict"

# Create the Variant Calls directory if it doesn't exist
mkdir -p $output_variant_dir

# Download the Singularity container
singularity pull $gatk_singularity_image

# Loop through recalibrated BAM files
for recal_bam in /sci/home/john.adeoye01/Desktop/mycoviruses/recalibrated/*/*.bam; do
    # Extract sample name from the recalibrated BAM file
    sample_name=$(basename "$(dirname "$recal_bam")")

    # Create a directory for each sample if it doesn't exist
    mkdir -p $output_variant_dir/$sample_name

    # Run GATK HaplotypeCaller
    singularity exec $gatk_singularity_image gatk --java-options "-Xmx4g" HaplotypeCaller \
        -I $recal_bam \
        -R $reference_genome \
        -O $output_variant_dir/$sample_name/${sample_name}_raw_variants.g.vcf.gz \
        -ploidy 1 \
        --sequence-dictionary $reference_dict

    # Optional: Index the gVCF file (useful for downstream analysis)
    singularity exec $gatk_singularity_image gatk IndexFeatureFile -I $output_variant_dir/$sample_name/${sample_name}_raw_variants.g.vcf.gz
done

# Combine variant calls across samples with GenotypeGVCFs
input_gvcfs=""
for gvcf in $output_variant_dir/*/*.g.vcf.gz; do
    input_gvcfs+="--variant $gvcf "
done

# Run GATK GenotypeGVCFs
singularity exec $gatk_singularity_image gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R $reference_genome \
    $input_gvcfs \
    -O $output_variant_dir/joint_variants.vcf.gz \
    --sequence-dictionary $reference_dict

# Index the joint VCF file
singularity exec $gatk_singularity_image gatk IndexFeatureFile -I $output_variant_dir/joint_variants.vcf.gz

