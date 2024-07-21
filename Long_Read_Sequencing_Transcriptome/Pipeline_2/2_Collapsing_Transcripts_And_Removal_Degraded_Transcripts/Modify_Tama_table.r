#!/usr/bin/Rscript

# R code for modifying table
# the transcript of the bed file needs to be replaced with the PB info
# Important to use the capped.bed file and not the trans_read.bed cause this one is transcripts before collapsing

library(tidyverse) 
filterdegradationtama<- read.delim("/Beatriz_Toledo/tama-master/filterdegradationtama.gtf", header=FALSE)
bff2nocap<-filterdegradationtama %>%
  mutate(V2=str_replace(filterdegradationtama$V2,"PBRI", "PacBio"))%>%
  mutate(V9=gsub('G', '"PB.', filterdegradationtama$V9)) 
bff2nocap<-bff2nocap %>%
  mutate(V9=str_replace(bff2nocap$V9,'exon_number ', 'exon_number "'))
bff2nocap<-bff2nocap %>%
  mutate(V9=gsub(';', '";', bff2nocap$V9)) write.table(bff2nocap,"/Beatriz_Toledo/tama-master/filterdegradationtama.modified.gtf", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
