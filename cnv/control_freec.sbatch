#!/bin/bash
# -----------------------------------------------------------------------------
# Description:        control-freec analysis
# Author:             Maria Laura Fabre (https://github.com/laurafabre)
# Date Created:       <2024-03-08>
# Last Modified:      <2024-03-08>
# Usage:              bash control_freec.sbatch 
# Notes:             It is needed to have a config file for each bam sample. Modify the paths especially in [sample] -> mateFile
# -----------------------------------------------------------------------------
#SBATCH --job-name=freec
#SBATCH --partition=standard 
#SBATCH --output=out/freec.out.%j
#SBATCH --error=out/freec.err.%j
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G 

# Load any necessary modules if needed
module purge

# activate env
mamba activate VCanalysis

# Define paths and filenames
#
# Get the directory path of the script
main_dir="freec/"
script_dir="freec/FREEC-11.6b/scripts/"

# Your commands to run Control-FREEC 
# Calculate the CNV

for subdir in "$main_dir"/*/; do
   # Check if the directory doesn't start with "FREEC", "reference", or "out"
   if [[ "$subdir" != FREEC* && "$subdir" != reference* && "$subdir" != out* ]]; then
       echo "Processing directory: $subdir"
       
       # Find the config file and check if it exists
       config_file=$(find "$subdir" -maxdepth 1 -type f -name "config_*.txt")
       if [ -f "$config_file" ]; then
           echo "Running FREEC for $config_file"
           FREEC-11.6b/src/freec -conf "$config_file"
       else
           echo "Config file not found in $subdir"
       fi
   fi
done


## Calculate the statistical significance of CNVs

# Loop through subdirectories
for subdir in "$main_dir"/*/; do
    echo "Processing directory: $subdir"
    if [ -d "$subdir" ]; then
        cd "$subdir" || exit 1
        
        # Run the R script for each pair of CNVs and ratio files
        for cnv_file in *_CNVs; do
            ratio_file="${cnv_file/_CNVs/_ratio.txt}"
            if [ -e "$cnv_file" ] && [ -e "$ratio_file" ]; then
                echo "Running R script for $cnv_file and $ratio_file"
                cat "$script_dir/assess_significance.R" | R --slave --args "$cnv_file" "$ratio_file"
            else
                echo "Files not found for $cnv_file and $ratio_file"
            fi
        done
    else
        echo "Directory $subdir does not exist or is inaccessible."
    fi
done
