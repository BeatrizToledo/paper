#!/usr/bin/env Rscript

#R code for identification of inclusion isoform of AS events

library(tidyr) 
library(tibble) 
library(plyr) 
library(dplyr) 
library(stringr) 
library(data.table)
library(mergeutils)
library(ggforce)
library(gridExtra)
library(ggplot2)
library(viridis)
library(ggbreak)
library(scales)
library(readxl)
library(unheadr)

keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
}

#filter out low expressed isoforms (less thant 4 reads in at least 1 splice-junction
jncalllistreadcount <- mergedpcstr2.nofus_NOCAP_SJsqanti_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
jncalllistreadcount1 <- jncalllistreadcount %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
jncalllistreadcountlessthan3 <- jncalllistreadcount1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
jncalllistreadcount1alljncmorethan3 <- anti_join(jncalllistreadcount1,jncalllistreadcountlessthan3)
inner_join(jncalllistreadcount1alljncmorethan3,jncalllistreadcountlessthan3)
transcriptsalljncmorethan3 <- jncalllistreadcount1alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3 <- mergedpcstr.nofus_NOCAP_sqanti_corrected %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=transcriptsalljncmorethan3, by.x = "V10", by.y = "isoform") 

#List with nodes
NODES <- rbind(PP.DP.down.merged,DP.N.down.merged) %>% rbind(DP.N.up.merged) %>% rbind(PP.DP.up.merged) %>% distinct() %>%
  separate(Coord, into=c("Chr","Coord"), sep=":") %>% separate(Coord, into=c("Start","End"), sep="-")

#CE
whipCEnotmerged <- NODES %>% merge(mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3,by.x = c("Start", "End", "Strand", "Chr"), by.y = c("V4", "V5", "V7","V1"))%>% subset(select = -c (V2, V3, V6, V8))
whipCEnotmergedname <- whipCEnotmerged%>%subset(select = c("name"))
all.notmerged.simplenamecheck <- NODES %>% subset(select = c("name")) %>% anti_join(whipCEnotmergedname) %>% merge(NODES, by.y="name")
#Maybe these are in a database
all.notmerged.simplenamecheck[grep('CE', all.notmerged.simplenamecheck$Type), ]



