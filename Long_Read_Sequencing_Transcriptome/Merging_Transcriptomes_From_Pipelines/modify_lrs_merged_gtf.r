#!/usr/bin/env Rscript

#R code for modifying .gtf
#should be in folder /../MERGED_TRANSCRIPTOME/
#the transcript G. needs to be replaced with the PB info

library(tidyverse)
library(dplyr) 

LRS_merged.gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_merged.gtf", header=FALSE,  quote="'")
LRS.gtf<-LRS_merged.gtf %>%
  mutate(V2=str_replace(LRS_merged.gtf$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('"G', '"PB.', LRS_merged.gtf$V9)) 
write.table(LRS.gtf,"/../MERGED_TRANSCRIPTOME/LRS.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
