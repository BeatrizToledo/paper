#!/usr/bin/Rscript

# R code for modifying table
# the transcript of the bed file needs to be replaced with the PB info
# Important to use the capped.bed file and not the trans_read.bed cause this one is transcripts before collapsing

library(stringr) 
library(dplyr) 
library(tidyverse)
library(plyr) 
 
filterdegradationtama<- read.delim("/../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.gtf", header=FALSE)
filterdegradationtama1 <-filterdegradationtama %>%
  mutate(V2=str_replace(filterdegradationtama$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('G', '"PB.', filterdegradationtama$V9)) 
filterdegradationtama1 <-filterdegradationtama1 %>%
  mutate(V9=str_replace(filterdegradationtama1$V9,'exon_number ', 'exon_number "'))
filterdegradationtama1<-filterdegradationtama1 %>%
  mutate(V9=gsub(';', '";', filterdegradationtama1$V9)) 
write.table(filterdegradationtama1,"/../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
