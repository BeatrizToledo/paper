#!/bin/bash

# Assembled transcripts might have artefacts or be low quality.
# Here some clean up to make before proceeding with SQANTI analysis.

# compare with mm10 and annotate the genes according to Ensembl

gffcompare /../ASSEMBLY/SRS.gtf -r mm10.filtered.novel.gtf -R -M -o /../GFFCOMPARE/SRS -V

# extract bad class code transcripts with awk from the annotated one
awk '{if($0~"class_code \"(s|e|x|p)\""){print $0}}' /../GFFCOMPARE/SRS.annotated.gtf > /../GFFCOMPARE/SRS.annotated.bad.tx.gtf

# run R script to remove bad quality transcripts and other
Rscript --vanilla << EOF

# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)
library(splitstackshape)

#read the GTF file with annotated bad ranscripts
SRS.annotated.bad.tx.gtf <- read.delim("/../GFFCOMPARE/SRS.annotated.bad.tx.gtf", header=FALSE)

gtf.bad <- SRS.annotated.bad.tx.gtf %>% 
  mutate(V9 = str_replace(V9, "\\;", "|")) %>% 
  separate(V9, into = c("transcript_id", "rest"), sep = "\\|")

# read the original GTF file
gtf <- read.delim("/../ASSEMBLY/SRS.gtf", header=FALSE)
gtf$ID <- paste(gtf$V9)

# split concatenated strings into multiple columns
gtf <- cSplit(gtf, "ID", sep= ";", type.convert=FALSE)

# subset the data to remove bad transcripts
gtf.good <- subset(gtf, !(gtf$ID_2 %in% gtf.bad$transcript_id))

# remove transcripts with only 1 bp length
unibp <- subset(gtf.good, V4 == V5)
gtf.good <- subset(gtf.good, !(gtf.good$ID_2 %in% unibp$ID_2))

# remove the transcripts with "."
gtf.good <- subset(gtf.good, V7 == "+" | V7 == "-")

# write the filtered data to a file
write.table(gtf.good, "/../ASSEMBLY/SRS.gtf", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)

EOF

# End of the shell script
