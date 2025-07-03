#!/bin/bash
#SBATCH --job-name=gatk_joint_variant_calling
#SBATCH --output=gatk_joint_variant_calling.out
#SBATCH --error=gatk_joint_variant_calling.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Load Singularity module
module load singularity

# Define paths
REFERENCE_GENOME="/sci/home/john.adeoye01/Desktop/mycoviruses/ref_genome/Afum293_genome.fasta"
OUTPUT_DIR="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/output"
ANC_GMM_VCF="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/anc_gmm/anc_gmm_raw_variants.g.vcf.gz"
RI_VCF="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/RI/RI_raw_variants.g.vcf.gz"
VC_CYCLOHEXAMIDE_VCF="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/VC-Cyclohexamide/VC-Cyclohexamide_raw_variants.g.vcf.gz"
VC_RIBAVIRIN_VCF="/sci/home/john.adeoye01/Desktop/mycoviruses/variant_calls/VC-Ribavirin/VC-Ribavirin_raw_variants.g.vcf.gz"

# Combine GVCFs
singularity run https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0 \
    gatk CombineGVCFs \
    -R ${REFERENCE_GENOME} \
    -V ${ANC_GMM_VCF} \
    -V ${RI_VCF} \
    -V ${VC_CYCLOHEXAMIDE_VCF} \
    -V ${VC_RIBAVIRIN_VCF} \
    -O ${OUTPUT_DIR}/combined_variants.g.vcf.gz

# Joint Genotyping
singularity exec https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0 \
    gatk GenotypeGVCFs \
    -R ${REFERENCE_GENOME} \
    -V ${OUTPUT_DIR}/combined_variants.g.vcf.gz \
    -O ${OUTPUT_DIR}/joint_variants.vcf.gz
