#!/bin/bash

# Assembled transcripts might have artefacts or be low quality.
# Here some clean up to make before proceeding with SQANTI analysis.

# Step 1: Compare with mm10 and annotate the genes according to Ensembl
gffcompare /../assembly/SRS.gtf -r SRS.gtf -R -M -o /../gffcompare/SRS -V

# Step 2: Extract bad class code transcripts with awk from the annotated one
awk '{if($0~"class_code \"(s|e|x|p)\""){print $0}}' /../gffcompare/SRS.annotated.gtf > /../gffcompare/SRS.annotated.bad.tx.gtf

# Step 3: Run R script to process the GTF files
Rscript --vanilla << EOF

# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)
library(splitstackshape)

# Read the bad transcripts GTF file
SRS.annotated.bad.tx.gtf <- read.delim("/../gffcompare/SRS.annotated.bad.tx.gtf", header=FALSE)

# Manipulate the data
gtf.bad <- SRS.annotated.bad.tx.gtf %>% 
  mutate(V9 = str_replace(V9, "\\;", "|")) %>% 
  separate(V9, into = c("transcript_id", "rest"), sep = "\\|")

# Read the original GTF file
gtf <- read.delim("/../assembly/SRS.gtf", header=FALSE)
gtf$ID <- paste(gtf$V9)

# Split concatenated strings into multiple columns
gtf <- cSplit(gtf, "ID", sep= ";", type.convert=FALSE)

# Subset the data to remove bad transcripts
gtf.good <- subset(gtf, !(gtf$ID_2 %in% gtf.bad$transcript_id))

# Remove transcripts with only 1 bp length
unibp <- subset(gtf.good, V4 == V5)
gtf.good <- subset(gtf.good, !(gtf.good$ID_2 %in% unibp$ID_2))

# Remove the transcripts with "."
gtf.good <- subset(gtf.good, V7 == "+" | V7 == "-")

# Write the filtered data to a file
write.table(gtf.good, "/../assembly/SRS.gtf", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)

EOF

# End of the shell script
