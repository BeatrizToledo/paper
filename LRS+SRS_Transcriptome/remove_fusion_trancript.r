#!/usr/bin/env Rscript

#R code for removing fusion trancripts
#should be in folder /../MERGED_TRANSCRIPTOME/

library(dplyr) 

LRS_SRS <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_classification.txt", header=FALSE,  quote="'")
LRS_SRS_nofusion_classification <- LRS_SRS %>% subset(structural_category != "fusion", select = c("isoform", "structural_category"))  
write.table(LRS_SRS_nofusion_classification,"/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion_classification.txt", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
