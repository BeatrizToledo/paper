#!/usr/bin/env Rscript

#R code for modifying .gtf
#should be in folder /../MERGED_TRANSCRIPTOME/
#the transcript G. needs to be replaced with the PB info

library(tidyverse)
library(dplyr) 

lrs_srs_merged.gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_merged.gtf", header=FALSE,  quote="'")
lrs_srs_merged_modified.gtf<-lrs_srs_merged.gtf %>%
  mutate(V2=str_replace(lrs_srs_merged.gtf$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('"G', '"PB.', lrs_srs_merged.gtf$V9)) 
write.table(lrs_srs_merged_modified.gtf,"/../MERGED_TRANSCRIPTOME/LRS_SRS_merged_modified.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
