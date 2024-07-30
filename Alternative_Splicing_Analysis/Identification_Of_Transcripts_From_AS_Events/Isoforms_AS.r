#!/usr/bin/env Rscript

#R code for identification of inclusion isoforms and exclusion isoforms of AS events

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

#filter out low expressed isoforms (less thant 4 reads in at least 1 splice-junction) from LRS_SRS

#table with splice junction coverage information
LRS_SRS_junctions <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion_junctions.txt", header=FALSE,  quote="'")
LRS_SRS_junctions_coverage <- LRS_SRS_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique")) %>% LRS_SRS_junctions_coverage %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)

#list of isoforms with at least 1 splice junction with less than 4 reads
LRS_SRS_junctions_coverage_lessthan4 <- LRS_SRS_junctions_coverage %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()

#list of isoforms with all splice junction with at least 4 reads
LRS_SRS_alljncmorethan3 <- anti_join(LRS_SRS_junctions_coverage1,LRS_SRS_junctions_coverage_lessthan4)
LRS_SRS_transcript_alljncmorethan3 <- LRS_SRS_alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

#gtf of isoforms with all splice junction with at least 4 reads
LRS_SRS_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion.gtf", header=FALSE,  quote="'")
LRS_SRS_gtf_alljncmorethan3 <- LRS_SRS_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=LRS_SRS_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 

#filter out low expressed isoforms (less thant 4 reads in at least 1 splice-junction) from ensembl and refseq

#table with splice junction coverage information
ensembl_junctions <- read.delim("/../REFERENCE_INPUT/ensembl_junctions.txt", header=FALSE,  quote="'")
refseq_junctions <- read.delim("/../REFERENCE_INPUT/refseq_junctions.txt", header=FALSE,  quote="'")
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

#gtf of isoforms with all splice junction with at least 4 reads
ensembl_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/ensembl_corrected.gtf", header=FALSE,  quote="'")
refseq_gtf <- read.delim("/../MERGED_TRANSCRIPTOME/refseq_corrected.gtf", header=FALSE,  quote="'")
ensembl_gtf_alljncmorethan3 <- ensembl_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=ensembl_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 
refseq_gtf_alljncmorethan3 <- refseq_gtf %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=refseq_transcript_alljncmorethan3, by.x = "V10", by.y = "isoform") 

#list with AS coordinates

NSC.NP.down <- read.delim("/../WHIPPET_RESULTS/NSC.NP.down.txt", header=FALSE,  quote="'")
NP.N.down <- read.delim("/../WHIPPET_RESULTS/NP.N.down.txt", header=FALSE,  quote="'")
NP.N.up <- read.delim("/../WHIPPET_RESULTS/NP.N.up.txt", header=FALSE,  quote="'")
NSC.NP.up <- read.delim("/../WHIPPET_RESULTS/NSC.NP.up.txt", header=FALSE,  quote="'")

#list with all AS coordinates
NODES <- rbind(NSC.NP.down,NP.N.down) %>% rbind(NSC.NP.up) %>% rbind(NSC.NP.up) %>% distinct() %>%
  separate(Coord, into=c("Chr","Coord"), sep=":") %>% separate(Coord, into=c("Start","End"), sep="-")

##########################

#check for inclusion isoforms from LRS_SRS transcriptome

#AS event is assigned to a inclusion isoform if their coordinates overlapped, at least partially, with an exon
#for this compare AS events coordinates with .gtf file
#AS events that did not have a corresponding inclusion isoform in LRS_SRS were checked for inclusion isoforms in ensembl and refseq

#check which AS events match with both Start and End coordinate of an exon
as_start_end <- NODES %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("Start", "End", "Strand", "Chr"), by.y = c("V4", "V5", "V7","V1"))%>% subset(select = -c (V2, V3, V6, V8)) %>% mutate (V9 = gsub('"', '', V9))
as_start_end_name <- as_start_end %>%subset(select = c("name"))
as_not_start_end <- NODES %>% subset(select = c("name")) %>% anti_join(as_start_end_name) %>% merge(NODES, by.y="name")

#check which of the remaining AS events that match with either Start or End coordinate of an exon
as_start <- as_not_start_end %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("Start", "Strand", "Chr"), by.y = c("V4", "V7","V1")) %>% subset(select = -c (V5, V2, V3, V6, V8)) %>% mutate (V9 = gsub('"', '', V9))
as_end <- as_not_start_end %>% merge(LRS_SRS_transcript_alljncmorethan3,by.x = c("End", "Strand", "Chr"), by.y = c("V5", "V7","V1"))%>% subset(select = -c (V4, V2, V3, V6, V8)) %>% mutate (V9 = gsub('"', '', V9))
as_start_name <- as_start %>%subset(select = c("name")) %>%distinct() 
as_end_name <-as_end%>%subset(select = c("name"))%>% distinct()
as_not_start_and_or_end <-  NODES %>% subset(select = c("name")) %>% setdiff(as_start_end_name) %>% setdiff(as_start_name) %>% setdiff(as_end_name) %>% merge(NODES, by="name")

#check which of the remaining AS events coordinates overlaps at least partially an exon
#prepare file for bedtools
as_not_start_and_or_end_coord <- as_not_start_and_or_end %>% mutate(Zero = 0) %>% subset(select =c("Chr", "Start", "End", "name", "Zero","Strand")) %>%  mutate(Chr=gsub('^', 'chr', Chr))
write.table(as_not_start_and_or_end_coord, "/../ISOFORM/as_not_start_and_or_end_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
LRS_SRS_alljncmorethan3_coord <- LRS_SRS_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(LRS_SRS_alljncmorethan3_coord, "/../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

system("/../ISOFORM/bedtools_inclusion_lrs_srs.sh")

#output from bedtools
as_overlap <- read.delim("/../ISOFORM/as_overlap.txt", header=FALSE,  quote="'")
as_overlap_name <- as_overlap %>%subset(select = c("V10")) %>% dplyr::rename(name=V10)
as_overlap <- as_overlap %>% subset(select = c("V1", "V4","V6", "V8", "V9", "V10")) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1", "V4", "V6", "V8", "V9", "V10"), new = c('Chr','V9', "Strand", "Start", "End", "name"))

#AS events without inclusion isoform in LRS_SRS
as_not_overlap <-  NODES %>% subset(select = c("name")) %>% distinct()%>% setdiff(as_start_end_name) %>% setdiff(as_start_name) %>% setdiff(as_end_name) %>%setdiff(as_overlap_name) %>%merge(NODES, by="name") 

#check if coordinate of the remaining AS events match coord from annotated transcripts

#check which of the remaining AS events match with both Start and End coordinate of an exon
ens_startend <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name")) %>% mutate (V9 = gsub('"', '', V9))
ens_startend_node <- ens_startend %>% subset(select = c("name")) %>% distinct()
refseq_startend <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name")) %>% anti_join(ens_startend_node) %>% mutate (V9 = gsub('"', '', V9))
refseq_startend_node <- refseq_startend %>% subset(select = c("name")) %>% distinct()

#check which of the remaining AS events have Start or End coordinates that match with either Start or End coordinate of an exon
ens_start <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(ens_startend_node) %>% anti_join(refseq_startend_node) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End) %>% mutate (V9 = gsub('"', '', V9))
ens_start_node <- ens_start %>% subset(select = c("name")) %>% distinct()
ens_end <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V5"), by.y = c("End"))  %>% anti_join(ens_startend_node) %>% anti_join(refseq_startend_node)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start) %>% mutate (V9 = gsub('"', '', V9))
ens_end_node <- ens_end %>% subset(select = c("name")) %>% distinct()
refseq_start <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V4"), by.y = c("Start")) %>% anti_join(ens_startend_node) %>% anti_join(refseq_startend_node) %>% anti_join(ens_start_node) %>% anti_join(ens_end_node) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End) %>% mutate (V9 = gsub('"', '', V9))
refseq_start_node <- refseq_start %>% subset(select = c("name")) %>% distinct()
refseq_end <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>% merge(as_not_overlap, by.x = c("V5"), by.y = c("End")) %>% anti_join(ens_startend_node) %>% anti_join(refseq_startend_node) %>% anti_join(ens_start_node) %>% anti_join(ens_end_node) %>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start) %>% mutate (V9 = gsub('"', '', V9))
refseq_end_node <- refseq_end %>% subset(select = c("name")) %>% distinct()

#check which of the remaining AS events overlaps an exon at least partially an exon
#prepare file for bedtools
as_not_overlap_coord <- as_not_overlap %>% mutate(Zero = 0) %>% subset(select =c("Chr", "Start", "End", "name", "Zero","Strand")) %>%  mutate(Chr=gsub('^', 'chr', Chr))
write.table(as_not_overlap_coord, "/../ISOFORM/as_not_overlap_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
ensembl_gtf_alljncmorethan3_coord <- ensembl_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(ensembl_gtf_alljncmorethan3_coord, "/../ISOFORM/ensembl_gtf_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
refseq_gtf_alljncmorethan3_coord <- refseq_gtf_alljncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(refseq_gtf_alljncmorethan3_coord, "/../ISOFORM/refseq_gtf_alljncmorethan3_coord.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

system("/../ISOFORM/bedtools_inclusion_ens_ref.sh")

as_overlap_ens <- read.delim("/../ISOFORM/as_overlap_ens.txt", header=FALSE)
as_overlap_ref <- read.delim("/../ISOFORM/as_overlap_ref.txt", header=FALSE)
as_overlap_ens_name <- as_overlap_ens %>%subset(select = c("V10")) %>% dplyr::rename(name=V10)
as_overlap_ens1 <- as_overlap_ens %>% subset(select = c("V1", "V4","V6", "V8", "V9", "V10")) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1", "V4", "V6", "V8", "V9", "V10"), new = c('Chr','V9', "Strand", "Start", "End", "name"))
as_overlap_ref_name <- as_overlap_ref %>%subset(select = c("V10")) %>% dplyr::rename(name=V10)
as_overlap_ref_name <- anti_join(as_overlap_ref_name,as_overlap_ens_name)
as_overlap_ref <- as_overlap_ref %>% subset(select = c("V1", "V4","V6", "V8", "V9", "V10")) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1", "V4", "V6", "V8", "V9", "V10"), new = c('Chr','V9', "Strand", "Start", "End", "name"))
as_overlap_ref1 <- anti_join(as_overlap_ref,as_overlap_ens_name)

#make one table table with the AS event and corresponding isoform
as_incl <- rbind(as_start_end, as_start) %>% rbind(as_end) %>% rbind(as_overlap) %>% rbind(ens_startend) %>% rbind(ens_start) %>% rbind(ens_end) %>% rbind(as_overlap_ens1) %>% rbind(refseq_startend) %>% rbind(refseq_start) %>% rbind(refseq_end) %>% rbind(as_overlap_ref1) 
as_inclusion_isoform <- as_incl %>% separate(V9, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr))  
as_inclusion_isoform_node <- as_inclusion_isoform %>% subset(select =c("name"))%>% distinct()

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
LRS_SRS_excl_event_tr <- read.delim("/../ISOFORM/LRS_SRS_excl_event_tr.txt", header=FALSE)
ensembl_excl_event_tr <- read.delim("/../ISOFORM/ensembl_excl_event_tr.txt", header=FALSE)
refseq_excl_event_tr <- read.delim("/../ISOFORM/refseq_excl_event_tr.txt", header=FALSE)

#exclusion isoforms in LRS_SRS
LRS_SRS_excl_nodes <- LRS_SRS_excl_event_tr %>% subset(select = c("V10")) %>% rename(name=V10) %>% distinct()
LRS_SRS_excl_nodes_iso <- LRS_SRS_excl_event_tr %>% subset(select = c("V4", "V10")) %>% distinct()

#AS events without exclusion isoform in LRS_SRS
excl_nodes_not_lrssrs <- anti_join(NODES,LRS_SRS_excl_nodes)

#check if coordinate of the remaining AS events is located within an intron of annotated transcripts
allensembl_excl_nodes <- ensembl_excl_event_tr %>% distinct() %>% subset(select = c("V10")) %>% rename(name=V10) %>% distinct()
allensembl_excl_nodes_iso <- ensembl_event_tr %>% distinct() %>% subset(select = c("V4", "V10")) %>% distinct()
ensembl_excl_nodes<- inner_join(excl_nodes_not_lrssrs,allensembl_excl_nodes)%>% distinct() %>% ensembl_excl_nodes %>% subset(select = c("name"))
ensembl_excl_nodes_iso <- ensembl_excl_nodes %>% distinct() %>% merge(allensembl_excl_nodes_iso, by.x = "name", by.y = "V10" )

#check in refseq
allrefseq_excl_nodes <- refseq_excl_event_tr %>% distinct() %>% subset(select = c("V10")) %>% rename(name=V10)
allrefseq_excl_nodes_iso <- refseq_event_tr %>% distinct() %>% subset(select = c("V4", "V10")) 
refseq_excl_nodes <- inner_join(excl_nodes_not_lrssrs,allrefseq_excl_nodes)%>% distinct() %>% subset(select = c("name"))
onlyrefseq_excl_nodes <- anti_join(refseq_excl_nodes,ensembl_excl_nodes)
refseq_excl_nodes_iso <- onlyrefseq_excl_nodes %>% distinct() %>% merge(allrefseq_excl_nodes_iso, by.x = "name", by.y = "V10" )

#join the exclusion isoform in one list
as_exclusion_isoform <- rbind(LRS_SRS_excl_nodes_iso,ensembl_excl_nodes_iso) %>% rbind(refseq_excl_nodes_iso) %>% distinct() %>% merge(NODES, by = "name")
as_exclusion_isoform_node <- as_exclusion_isoform %>% subset(select =c("name")) %>% distinct()

####################################

# select nodes assigned to at least one inclusion and one exclusion isoform for further analysis
NODESTOBEANALYZED <- intersect(as_exclusion_isoform_node,as_inclusion_isoform_node) %>% distinct()

as_exclusion_isoform_tobeanalyzed <- as_exclusion_isoform %>% merge(NODESTOBEANALYZED, by = "name")
as_inclusion_isoform_tobeanalyzed <- as_inclusion_isoform %>% merge(NODESTOBEANALYZED, by = "name")




