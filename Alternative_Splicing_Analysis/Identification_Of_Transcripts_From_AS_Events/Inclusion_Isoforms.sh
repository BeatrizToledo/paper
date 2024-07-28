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

#filter out low expressed isoforms (less thant 4 reads in at least 1 splice-junction)

#table with splice junction coverage information
LRS_SRS_junctions <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion_junctions.txt", header=FALSE,  quote="'")
LRS_SRS_junctions_coverage <- LRS_SRS_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
LRS_SRS_junctions_coverage1 <- LRS_SRS_junctions_coverage %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)

#list of isoforms with at least 1 splice junction with less than 4 reads
LRS_SRS_junctions_coverage_lessthan4 <- LRS_SRS_junctions_coverage1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()

#list of isoforms with all splice junction with at least 4 reads
LRS_SRS_alljncmorethan3 <- anti_join(LRS_SRS_junctions_coverage1,LRS_SRS_junctions_coverage_lessthan4)
LRS_SRS_transcript_alljncmorethan3 <- LRS_SRS_alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

LRS_SRS_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion.gtf", header=FALSE,  quote="'")
LRS_SRS_gtf_alljncmorethan3 <- LRS_SRS_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=LRS_SRS_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 

#list with AS coordinates
NSC.NP.down <- read.delim("/../WHIPPET_RESULTS/NSC.NP.down.txt", header=FALSE,  quote="'")
NP.N.down <- read.delim("/../WHIPPET_RESULTS/NP.N.down.txt", header=FALSE,  quote="'")
NP.N.up <- read.delim("/../WHIPPET_RESULTS/NP.N.up.txt", header=FALSE,  quote="'")
NSC.NP.up <- read.delim("/../WHIPPET_RESULTS/NSC.NP.up.txt", header=FALSE,  quote="'")

#list with all AS coordinates
NODES <- rbind(NSC.NP.down,NP.N.down) %>% rbind(NSC.NP.up) %>% rbind(NSC.NP.up) %>% distinct() %>%
  separate(Coord, into=c("Chr","Coord"), sep=":") %>% separate(Coord, into=c("Start","End"), sep="-")

##########################

#Check for inclusion isoforms from LRS_SRS transcriptome

#Check which AS events have Start and End coordinates that match with both Start and End coordinate of an exon
as_start_end <- NODES %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("Start", "End", "Strand", "Chr"), by.y = c("V4", "V5", "V7","V1"))%>% subset(select = -c (V2, V3, V6, V8))
as_start_end_name <- as_start_end %>%subset(select = c("name"))
as_not_start_end <- NODES %>% subset(select = c("name")) %>% anti_join(as_start_end_name) %>% merge(NODES, by.y="name")

#Check which of the remaining AS events have Start or End coordinates that match with either Start or End coordinate of an exon
as_start <- as_not_start_end %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("Start", "Strand", "Chr"), by.y = c("V4", "V7","V1")) %>% subset(select = -c (V5, V2, V3, V6, V8)) 
as_end <- as_not_start_end %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("End", "Strand", "Chr"), by.y = c("V5", "V7","V1"))%>% subset(select = -c (V4, V2, V3, V6, V8)) 
as_start_name <- as_start %>%subset(select = c("name")) %>%distinct() 
as_end_name <-as_end%>%subset(select = c("name"))%>% distinct()
as_not_start_and_or_end <-  NODES %>% subset(select = c("name")) %>% setdiff(as_start_end_name) %>% setdiff(as_start_name) %>% setdiff(as_end_name) %>% merge(NODES, by="name")

#check which of the remaining AS events coordinates overlaps an exon
#prepare file for bedtools
as_not_start_and_or_end1 <- as_not_start_and_or_end %>% mutate(Zero = 0) %>% subset(select =c("Chr", "Start", "End", "name", "Zero","Strand")) %>%  mutate(Chr=gsub('^', 'chr', Chr))
write.table(as_not_start_and_or_end1, "/../ISOFORM/as_not_start_and_or_end_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

LRS_SRS_alljncmorethan3_coord <- LRS_SRS_gtf_alljncmorethan3 %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>% distinct() %>%
  separate(V9, into = c("a", "b", "c"), sep = ";"  ) %>% unite(V9, c("b", "a", "c"), sep = "; ") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>%
  subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(LRS_SRS_alljncmorethan3_coord, "/../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

#on shell do bedtools intersect -a /../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt -b /../ISOFORM/as_not_start_and_or_end_coord.txt -wo -s -F 0.99 -filenames > /../ISOFORM/as_overlap.txt
system("/../ISOFORM/bedtools_lrs_srs.sh")

as_overlap <- read.delim("/../ISOFORM/as_overlap.txt", header=FALSE,  quote="'")
as_overlap_name <- as_overlap %>%subset(select = c("V10")) %>% dplyr::rename(name=V10)
as_overlap <- as_overlap %>% subset(select = c("V1", "V4","V6", "V8", "V9", "V10")) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1", "V4", "V6", "V8", "V9", "V10"), new = c('Chr','V9', "Strand", "Start", "End", "name"))

#########################





#########################

EVENTS <- whipCEnotmerged%>% rbind(whipAAADnotmerged) %>%rbind(whipAAAD1notmerged) %>%subset( select= c('Chr','V9', "Strand", "Start", "End", "name")) %>%rbind(notmergedoverlatirwhip1) %>% distinct() %>% mutate (V9 = gsub('"', '', V9))
NODES %>% subset(select= c("name")) %>% distinct() %>%nrow()
EVENTSname <- EVENTS %>% subset(select= c("name")) %>% distinct()
#total 2682, 68 missing

#AS events that did not overlap 
as_not_overlap <-  NODES %>% subset(select = c("name")) %>% distinct()%>% setdiff(as_start_end_name) %>% setdiff(as_start_name) %>% setdiff(as_end_name) %>%setdiff(as_overlap_name) %>%merge(NODES, by="name") 





#check if coordinate of the remaining AS events match coord from annotated transcripts
#filter out low expressed isoforms (less thant 4 reads in at least 1 splice-junction)

#table with splice junction coverage information
ensembl_junctions <- read.delim("/../REFERENCE_INPUT/ensembl_junctions.txt", header=FALSE,  quote="'")
refseq_junctions <- read.delim("/../REFERENCE_INPUT/refseq_junctions.txt", header=FALSE,  quote="'")
ensembl_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/ensembl_corrected.gtf", header=FALSE,  quote="'")
refseq_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/refseq_corrected.gtf", header=FALSE,  quote="'")

ensembl_junctions_coverage <- ensembl_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique")) %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
refseq_junctions_coverage <- refseq_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique")) %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)

#list of isoforms with at least 1 splice junction with less than 4 reads
ensembl_junctions_coverage_lessthan4 <- ensembl_junctions_coverage1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
refseq_junctions_coverage_lessthan4 <- refseq_junctions_coverage1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()

#list of isoforms with all splice junction with at least 4 reads
ensembl_alljncmorethan3 <- anti_join(ensembl_junctions_coverage,ensembl_junctions_coverage_lessthan4)
refseq_alljncmorethan3 <- anti_join(refseq_junctions_coverage,refseq_junctions_coverage_lessthan4)
ensembl_transcript_alljncmorethan3 <- ensembl_alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()
refseq_transcript_alljncmorethan3 <- refseq_alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

#gtf (exon coordinates)
ensembl_gtf_alljncmorethan3 <- ensembl_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=ensembl_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 
refseq_gtf_alljncmorethan3 <- refseq_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=refseq_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 

#Check which of the remaining AS events have Start and End coordinates that match with both Start and End coordinate of an exon
ens_startend <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name"))
#PB.9743-19, PB.9582-18,PB.9316-9, PB.8772-6,PB.2280-4,  PB.12316-82, PB.4866-14 - all ens
refseq_startend <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name"))
refseq_startend <- refseq_startend[grep('PB.4732-34|PB.1991-12', refseq_startendE$name), ]
#PB.4732-34,PB.1991-12

#Check which of the remaining AS events have Start or End coordinates that match with either Start or End coordinate of an exon
ens_start <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(ens_startend) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End)
#PB.6167-3, PB.5031-60, PB.4800-6, PB.4866-7, PB.2781-9, PB.8482-32 - all ens
ens_end <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V5"), by.y = c("End"))  %>% anti_join(ens_startend)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start)
#PB.10632-10, PB.2400-2, PB.10531-13 
ens_end <- ens_end[-grep('PB.4732-34', ens_end$name), ]
refseq_start <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(ens_startend) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End)
refseq_start <- refseq_start[grep('PB.4052-20|PB.5326-30|PB.12524-23|PB.3724-21|PB.311-28|PB.1400-48|PB.449-3', refseq_start$name), ]
#PB.4052-20|PB.5326-30|PB.12524-23|PB.3724-21|PB.311-28|PB.1400-48|PB.449-3
refseq_end <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V5"), by.y = c("End"))  %>% anti_join(ens_startend)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start)
#PB.10942-29|PB.8463-5|PB.10400-12
refseq_end <- refseq_end[grep('PB.10942-29|PB.8463-5|PB.10400-12', refseq_end$name), ]

#check which of the remaining AS events coordinates overlaps an exon
#prepare file for bedtools
ensembl_gtf_alljncmorethan3_coord <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(ensembl_gtf_alljncmorethan3_coord, "/Beatriz_Toledo/nocapMerge/ensembl_gtf_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
refseq_gtf_alljncmorethan3_coord <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(refseq_gtf_alljncmorethan3_coord, "/Beatriz_Toledo/nocapMerge/refseq_gtf_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

#on shell do bedtools intersect -a LRS_SRS_alljncmorethan3_coord.txt -b as_not_start_and_or_end_coord.txt -wo -s -F 0.99 -filenames > as_overlap.txt
system("/../__/custom_script.sh")
#bedtools intersect -a ensembl_gtf_alljncmorethan3_coord.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > NODES2022.merged.overlatirensemblFilteredincluded.txt
#bedtools intersect -a refseq_gtf_alljncmorethan3_coord.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > NODES2022.merged.overlatirrefseqFilteredincluded.txt

notmergedoverlatirensembl <- read.delim("/Beatriz_Toledo/nocapMerge/NODES2022.merged.overlatirensemblFilteredincluded.txt", header=FALSE)
notmergedoverlatirrefseq <- read.delim("/Beatriz_Toledo/nocapMerge/NODES2022.merged.overlatirrefseqFilteredincluded.txt", header=FALSE)
reftocheck <- merge(all.notmerged.simplenamecheck4, notmergedoverlatirrefseq, by.x = "name", by.y = "V10")
reftocheck1<- reftocheck[grep('PB.11846-27|PB.4956-74', reftocheck$name), ] %>% subset(select = c("V1", "Start", "End", "Strand", "V4", "name"))
#PB.11846-27|PB.4956-74
enstocheck <- merge(all.notmerged.simplenamecheck4, notmergedoverlatirensembl, by.x = "name", by.y = "V10") 
#PB.13554-20|PB.13554-21,22|PB.1837-3|PB.2933-14|PB.2952-10|PB.4067-54|PB.9574-12
enstocheck1<- enstocheck[grep('PB.13554-20|PB.13554-21,22|PB.1837-3|PB.2933-14|PB.2952-10|PB.4067-54|PB.9574-12', enstocheck$name), ]%>% subset(select = c("V1", "Start", "End", "Strand", "V4", "name"))
ANOThave <- rbind(all.ensCE,all.ensstart) %>% rbind(all.ensend1) %>% rbind(all.RefSeqCE1) %>% rbind(all.RefSeqstart1) %>% rbind(all.RefSeqend1)  %>% subset(select =c('name'))
ANOThave1 <- rbind(reftocheck1) %>% rbind(enstocheck1) %>% subset(select =c('name'))
ANOThave2 <- rbind(all.ensCE,all.ensstart) %>% rbind(all.ensend1) %>% rbind(all.RefSeqCE1) %>% rbind(all.RefSeqstart1) %>% rbind(all.RefSeqend1)  
ANOThave3 <- rbind(reftocheck1) %>% rbind(enstocheck1) 
TOCHECK <-anti_join(all.notmerged.simplenamecheck4, ANOThave ) %>% anti_join( ANOThave1 ) 
#PB.12450-7, PB.13189-16, PB.235-73 - THESE 3 are in fact new
#28 not have a corresponding transcript with reasonable sj levels - 31 2682

EVENTSx<-EVENTS[grep('; tr', EVENTS$V9), ] %>% separate(V9, into = c("gene", "tr", "x"), sep = ";") %>% mutate (tr = gsub(' transcript_id ', '', tr))
EVENTSx2<-EVENTS[-grep('; tr', EVENTS$V9), ] %>% separate(V9, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr))
ANOThave2x2<-ANOThave2 %>% separate(V9, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr))  %>% setnames( old = c("V1", "V7", "V4", "V5"), new = c("Chr","Strand", "Start", "End"))
ANOThave3x2<-ANOThave3 %>% separate(V4, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr)) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1"), new = c('Chr'))
EVENTSX <- rbind(EVENTSx,EVENTSx2)%>% rbind(ANOThave2x2) %>% rbind(ANOThave3x2)%>% mutate (tr = gsub('"', '', tr))
EVENTSXnodes <- EVENTSX %>% subset(select =c("name"))%>% distinct()
NODESTOBEANALYZED <- intersect(nodesisoexcnodes,EVENTSXnodes)
#2581 NODES FOR NEXT ANALYSES, But have to recheck removing the fusion opnes


cdscoord <- mergedpcstr2.nofus_NOCAP_sqanti_corrected.gtf.cds %>% subset(V3 == "CDS") %>% separate(V9, into = c("V9", "V10"), sep = ";") %>% mutate (V10 = gsub(' transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10))
cdscoord1 <- cdscoord %>% dplyr::group_by(V10) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V9", "V10","V7", "start", "end")) %>% distinct()
cdscoordens <- Ensembl_corrected.gtf.cds %>% subset(V3 == "CDS") %>% separate(V9, into = c("V9", "V10", "V11"), sep = ";") %>% mutate (V10 = gsub(' transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10))
cdscoordens1 <- cdscoordens %>% dplyr::group_by(V10) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V9", "V10","V7", "start", "end")) %>% distinct()
cdscoordref <- mm10.ncbiRefSeq2 %>% subset(V3 == "CDS")  %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V11 = gsub('transcript_id ', '', V11)) %>% mutate (V11 = gsub('"', '', V11)) 
cdscoordref1 <- cdscoordref %>% dplyr::group_by(V11) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V10", "V11","V7", "start", "end")) %>% distinct() %>% rename(V9 = V10) %>% rename(V10 = V11)
cdscoordall<- cdscoord1 %>% rbind(cdscoordens1) %>% rbind(cdscoordref1)


################################

#check for exclusion isoforms from LRS_SRS transcriptome

#AS event is assigned to a exclusion isoform if it is located within an intron
#for this check if the AS events coordinates is located between the "genomic_start_coord" and "genomic_end_coord" from junction file with bedtools
#AS events that did not have a corresponding exclusion isoform in LRS_SRS were checked for exclusion isoforms in ensembl and refseq

#prepare file for bedtools
NODES_coord <- NODES%>%mutate(Zero = 0) %>%  mutate(Chr=gsub('^', 'chr', Chr))%>% subset(select =c("Chr", "Start", "End","name", "Zero" ,"Strand")) %>% subset(End>Start) %>% distinct()
ensembl_alljncmorethan3_coord <- ensembl_alljncmorethan3 %>% mutate(Zero = 0) %>% subset(select =c("chrom", "genomic_start_coord", "genomic_end_coord", "isoform", "Zero","strand")) %>%  mutate(chrom=gsub('^', 'chr', chrom)) 
refseq_alljncmorethan3_coord <- refseq_alljncmorethan3 %>% mutate(Zero = 0) %>% subset(select =c("chrom", "genomic_start_coord", "genomic_end_coord", "isoform", "Zero","strand")) %>%  mutate(chrom=gsub('^', 'chr', chrom)) 
LRS_SRS_alljncmorethan3_coord <- LRS_SRS_alljncmorethan3 %>% mutate(Zero = 0) %>% subset(select =c("chrom", "genomic_start_coord", "genomic_end_coord", "isoform", "Zero","strand")) %>%  mutate(chrom=gsub('^', 'chr', chrom)) 

write.table(ensembl_alljncmorethan3_coord, "/../ISOFORM/ensembl_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
write.table(refseq_alljncmorethan3_coord, "/../ISOFORM/refseq_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
write.table(LRS_SRS_alljncmorethan3_coord, "/../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
write.table(NODES_coord, "/../ISOFORM/NODES_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

 system("/../ISOFORM/bedtools_exclusion.sh")

#output from bedtools
eventsexctrpcstr.2022.merged <- read.delim("LRS_SRS_excl_event_tr.txt", header=FALSE)
eventsexctrens.2022.merged <- read.delim("/ensembl_excl_event_tr.txt", header=FALSE)
eventsexctrrefseq.2022.merged <- read.delim("refseq_excl_event_tr.2022.merged.txt", header=FALSE)





