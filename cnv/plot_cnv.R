# -----------------------------------------------------------------------------
# Description:        plot supplementary figure
# Author:             Maria Laura Fabre (https://github.com/laurafabre)
# Date Created:       <2024-10-01>
# Last Modified:      <2024-10-01>
# -----------------------------------------------------------------------------
library(karyoploteR)
library(GenomicRanges)
library(dplyr)
library(CopyNumberPlots)

# In Control-FREEC’s *_ratio.txt file, a Ratio value of -1 indicates that the ratio could not be calculated for that particular position, usually due to insufficient coverage or other quality control issues.
data <- read.table('VC_Riba_freec2circos.txt')
#data <- read.table('VC_Cycloh_freec2circos.txt')
# Filter out rows where the value in column V4 is -1
data <- subset(data, V4 != -1)

# Calculate the log2 of the filtered data with an offset
data$lRR <- log2(data$V4 + 0.01)
data$V1 <- sub("Chr([0-9]+)_.*", "chr\\1", data$V1)
data$V1 <- sub("hs", "", data$V1)
#data$V4 <- NULL
colnames(data) <- c("Chromosome", 'Start', "End", "MedianRatio", "lrr")

# Load your karyotype file
karyotype_file <- read.table("supplementary plot/karyotype.txt"  )# Replace with the path to your karyotype file
colnames(karyotype_file) <- c("Chromosome", 'name', 'start', "end", "length")


# Convert to GRanges object
tiles <- makeGRangesFromDataFrame(data, 
                                  keep.extra.columns = TRUE, 
                                  seqnames.field = "Chromosome",
                                  start.field = "Start",
                                  end.field = "End",
                                  ignore.strand = TRUE)

# Check the GRanges object
tiles
tiles$lrr <- data$lrr

# Create GRanges object for the copy number regions
rr <- GRanges(
  seqnames = data$Chromosome,
  ranges = IRanges(start = data$Start, end = data$End),
  cn = data$MedianRatio,
  lrr= data$lrr# Use 'lrr' or another appropriate column for copy number calls
)

rr

# Define custom chromosomes and their lengths
chromosomes <- data %>%
  group_by(Chromosome) %>%
  summarize(length = max(End))


library(ggplot2)
library(dplyr)

plotStackedCopyNumbersGG <- function(data_files, karyotype_file, plot_title, colors) {
  
  all_data <- list()
  
  # Step 1: Load and preprocess each dataset
  for (i in seq_along(data_files)) {
    data <- read.table(data_files[i])  # Fix: Use data_files[i] to read the correct file
    
    # Filter out rows where the value in column V4 is -1
    data <- subset(data, V4 != -1)
    
    # Calculate the log2 of the filtered data with an offset
    data$lRR <- log2(data$V4 + 0.01)
    #as there are no significant CNV in common after running controlFREEC and DELLy we define cnv=1
    # Replace the cn column with 1 for all rows
    #data$cn <- data$V4
    data$cn <- 1
    data$V1 <- sub("Chr([0-9]+)_.*", "chr\\1", data$V1)
    data$V1 <- sub("hs", "", data$V1)
    colnames(data) <- c("Chromosome", 'Start', "End", "MedianRatio", "lrr", "cn")
    
    # Add a column to indicate which dataset this is
    #data$Dataset <- paste("Dataset", i)
    # Extract the basename of the file
    file_name <- basename(data_files[i])
    
    # Set the Dataset name based on the file name pattern
    if (grepl("RI", file_name)) {
      dataset_name <- "RI"
    } else if (grepl("VC_Cycloh", file_name)) {
      dataset_name <- "VC_Cyclohexamide"
    } else if (grepl("VC_Riba", file_name)) {
      dataset_name <- "VC_Ribavirin"
    } else {
      dataset_name <- paste("Dataset", i)  # Default name if no match
    }
    
    # Assign the real name to the Dataset column
    data$Dataset <- dataset_name
    
    all_data[[i]] <- data
  }
  
  # Combine all datasets into one data frame
  combined_data <- do.call(rbind, all_data)
  # Ensure the levels of the 'Dataset' factor have 'RI' as the last level
  combined_data$Dataset <- factor(combined_data$Dataset, levels = c("VC_Ribavirin", "VC_Cyclohexamide", "RI"))
  
  # Step 2: Load the karyotype data
  karyo_data <- read.table(karyotype_file, header = TRUE)
  colnames(karyo_data) <- c("Chromosome", 'name', 'start', "end", "length")
  
  # Convert lengths from base pairs to megabases
  karyo_data$length_MB <- karyo_data$length / 1e6
  
  # Step 3: Prepare karyotype positions for plotting
  combined_data <- combined_data %>%
    left_join(karyo_data, by = "Chromosome") %>%
    mutate(Start_MB = Start / 1e6,   # Convert Start to Mb
           End_MB = End / 1e6)       # Convert End to Mb
  
  # Step 4: Create the ggplot with faceting by Chromosome and Dataset
  p <- ggplot(combined_data, aes(x = Start_MB, y = cn, group = Dataset, color = Dataset)) +
    geom_line() +
    facet_grid(Dataset ~ Chromosome, scales = "free_x", space = "free_x") +
    scale_color_manual(values = colors) +
    scale_y_continuous(name = " Copy number", breaks = seq(0, 3, by = 1), limits = c(0, 2)) +  # Fix y-axis breaks
    scale_x_continuous(name = "Genomic Position (Mb)", labels = scales::comma_format()) +  # Update x-axis label
    theme_bw() +
    theme(
      text = element_text(size = 16),  # Sets the base size for all text
      axis.text = element_text(size = 14))+
    labs(title = "", y = "") 
  
  # Print the plot
  print(p)
}

# Example usage with 3 files and 3 colors:
plotStackedCopyNumbersGG(
  data_files = c( 'RI_freec2circos.txt',
                  'VC_Riba_freec2circos.txt',
                 'VC_Cycloh_freec2circos.txt'
                ),
  karyotype_file = 'supplementary plot/karyotype.txt',
  plot_title = "Stacked Fungal Genome Data",
  colors = c("darkblue", "orange", "darkred")
)

# 


