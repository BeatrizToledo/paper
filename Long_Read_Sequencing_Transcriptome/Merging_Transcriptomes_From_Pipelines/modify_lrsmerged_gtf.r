#!/usr/bin/env Rscript

#R code for modifying .gtf
#should be in folder /../MERGED_TRANSCRIPTOME/
#the transcript G. needs to be replaced with the PB info

library(tidyverse)
library(dplyr) 

LRSmerged.gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRSmerged.gtf", header=FALSE,  quote="'")
LRSmerged_modified.gtf<-LRSmerged.gtf %>%
  mutate(V2=str_replace(LRSmerged.gtf$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('"G', '"PB.', LRSmerged.gtf$V9)) 
write.table(LRSmerged_modified.gtf,"/../MERGED_TRANSCRIPTOME/LRSmerged_modified.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
