#!/usr/bin/env Rscript

#R code for modifying .gtf
#should be in folder /../MERGED_TRANSCRIPTOME/
#the transcript G. needs to be replaced with the PB info

library(tidyverse)
library(dplyr) 

LRS_SRS.gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS.gtf", header=FALSE,  quote="'")
LRS_SRS.gtf<-LRS_SRS.gtf %>%
  mutate(V2=str_replace(LRS_SRS.gtf$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('"G', '"PB.', LRS_SRS.gtf$V9)) 
write.table(LRS_SRS.gtf,"/../MERGED_TRANSCRIPTOME/LRS_SRS.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
