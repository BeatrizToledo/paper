#!/usr/bin/env Rscript

#R code for modifying table
#should be in folder /../TAMA/
#the transcript of the bed file needs to be replaced with the PB info
#important to use the nocap.bed file and not the trans_read.bed. trans_read.bed are the transcripts before collapsing

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
