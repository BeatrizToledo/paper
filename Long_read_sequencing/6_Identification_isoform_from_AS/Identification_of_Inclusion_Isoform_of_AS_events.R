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

#First step is to filterthe cases that have less than 4 sj

jncalllistreadcount <- mergedpcstr2.nofus_NOCAP_SJsqanti_junctions %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
jncalllistreadcount1 <- jncalllistreadcount %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
jncalllistreadcountlessthan3 <- jncalllistreadcount1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
jncalllistreadcount1alljncmorethan3 <- anti_join(jncalllistreadcount1,jncalllistreadcountlessthan3)
inner_join(jncalllistreadcount1alljncmorethan3,jncalllistreadcountlessthan3)
transcriptsalljncmorethan3 <- jncalllistreadcount1alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

ncbiRefSeq2jncalllistreadcount <- ncbiRefSeq2_SJsqanti %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
ncbiRefSeq2jncalllistreadcount1 <- ncbiRefSeq2jncalllistreadcount %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
ncbiRefSeq2jncalllistreadcountlessthan3 <- ncbiRefSeq2jncalllistreadcount1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
ncbiRefSeq2jncalllistreadcount1alljncmorethan3 <- anti_join(ncbiRefSeq2jncalllistreadcount1,ncbiRefSeq2jncalllistreadcountlessthan3)
inner_join(ncbiRefSeq2jncalllistreadcount1alljncmorethan3,ncbiRefSeq2jncalllistreadcountlessthan3)
ncbiRefSeq2transcriptsalljncmorethan3 <- ncbiRefSeq2jncalllistreadcount1alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

gencodejncalllistreadcount <- gencode_SJsqanti %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
gencodejncalllistreadcount1 <- gencodejncalllistreadcount %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
gencodejncalllistreadcountlessthan3 <- gencodejncalllistreadcount1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
gencodejncalllistreadcount1alljncmorethan3 <- anti_join(gencodejncalllistreadcount1,gencodejncalllistreadcountlessthan3)
inner_join(gencodejncalllistreadcount1alljncmorethan3,gencodejncalllistreadcountlessthan3)
gencodetranscriptsalljncmorethan3 <- gencodejncalllistreadcount1alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

ensembljncalllistreadcount <- ensembl_SJsqanti %>% subset(select = c("isoform", "genomic_start_coord", "genomic_end_coord", "chrom", "strand", "total_coverage_unique"))
ensembljncalllistreadcount1 <- ensembljncalllistreadcount %>%mutate(genomic_start_coord = genomic_start_coord-1) %>%   mutate(genomic_end_coord = genomic_end_coord+1)
ensembljncalllistreadcountlessthan3 <- ensembljncalllistreadcount1 %>% subset(total_coverage_unique <= 3) %>% subset(select = c("isoform")) %>% distinct()
ensembljncalllistreadcount1alljncmorethan3 <- anti_join(ensembljncalllistreadcount1,ensembljncalllistreadcountlessthan3)
inner_join(ensembljncalllistreadcount1alljncmorethan3,ensembljncalllistreadcountlessthan3)
ensembltranscriptsalljncmorethan3 <- ensembljncalllistreadcount1alljncmorethan3  %>% subset(select = c("isoform")) %>% distinct()

ensembl_SJsqanti_correctedjncmorethan3 <- ensembl_SJsqanti_corrected %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=ensembltranscriptsalljncmorethan3, by.x = "V10", by.y = "isoform") 
gencode_SJsqanti_correctedjncmorethan3 <- gencode_SJsqanti_corrected %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=gencodetranscriptsalljncmorethan3, by.x = "V10", by.y = "isoform") 
ncbiRefSeq2_SJsqanti_correctedjncmorethan3 <- ncbiRefSeq2_SJsqanti_corrected %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=ncbiRefSeq2transcriptsalljncmorethan3, by.x = "V10", by.y = "isoform") 
mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3 <- mergedpcstr.nofus_NOCAP_sqanti_corrected %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>%
  keeping.order(merge, y=transcriptsalljncmorethan3, by.x = "V10", by.y = "isoform") 

PP.DP.down.2022.merged %>% distinct() %>% nrow()
#262
PP.DP.up.2022.merged %>% distinct() %>% nrow()
#423
DP.N.down.2022.merged %>% distinct() %>% nrow()
#916
DP.N.up.2022.merged %>% distinct() %>% nrow()
#1417

PP.DP.down.2022.merged %>% merge(DP.N.up.2022.merged, by = "name")
#17
PP.DP.down.2022.merged %>% merge(DP.N.down.2022.merged, by = "name")
#75
PP.DP.up.2022.merged %>% merge(DP.N.up.2022.merged, by = "name")
#227
PP.DP.up.2022.merged %>% merge(DP.N.down.2022.merged, by = "name")
#17

NODES <- rbind(PP.DP.down.2022.merged,DP.N.down.2022.merged) %>% rbind(DP.N.up.2022.merged) %>% rbind(PP.DP.up.2022.merged) %>% distinct() %>%
  separate(Coord, into=c("Chr","Coord"), sep=":") %>% separate(Coord, into=c("Start","End"), sep="-")

#CE
whipCEnotmerged <- NODES %>% merge(mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3,by.x = c("Start", "End", "Strand", "Chr"), by.y = c("V4", "V5", "V7","V1"))%>% subset(select = -c (V2, V3, V6, V8))
whipCEnotmergedname <- whipCEnotmerged%>%subset(select = c("name"))
all.notmerged.simplenamecheck <- NODES %>% subset(select = c("name")) %>% anti_join(whipCEnotmergedname) %>% merge(NODES, by.y="name")
#Maybe these are in a database
all.notmerged.simplenamecheck[grep('CE', all.notmerged.simplenamecheck$Type), ]

#AA AD
whipAAADnotmerged <- all.notmerged.simplenamecheck %>% merge(mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3,by.x = c("Start", "Strand", "Chr"), by.y = c("V4", "V7","V1")) %>% subset(select = -c (V5, V2, V3, V6, V8)) 
whipAAAD1notmerged <- all.notmerged.simplenamecheck %>% merge(mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3,by.x = c("End", "Strand", "Chr"), by.y = c("V5", "V7","V1"))%>% subset(select = -c (V4, V2, V3, V6, V8)) 
whipAAADnotmergedname <- whipAAADnotmerged%>%subset(select = c("name")) %>%distinct() 
whipAAAD1notmergedname <-whipAAAD1notmerged%>%subset(select = c("name"))%>% distinct()
all.notmerged.simplenamecheck2<-  NODES %>% subset(select = c("name")) %>% setdiff(whipCEnotmergedname) %>% setdiff(whipAAADnotmergedname) %>% setdiff(whipAAAD1notmergedname) %>% merge(NODES, by="name")
all.notmerged.simplenamecheck2[grep('AA', all.notmerged.simplenamecheck2$Type), ]

#IR
all.notmerged.simplenamecheck3 <- all.notmerged.simplenamecheck2 %>% mutate(Zero = 0) %>% subset(select =c("Chr", "Start", "End", "name", "Zero","Strand")) %>%  mutate(Chr=gsub('^', 'chr', Chr))
write.table(all.notmerged.simplenamecheck3, "/Beatriz_Toledo/nocapMerge/NODES2022.merged.eventstocheckIR.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
gtfall <- mergedpcstr.nofus_NOCAP_sqanti_correctedjncmorethan3 %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V10 = gsub('transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10)) %>% distinct() %>%
  separate(V9, into = c("a", "b", "c"), sep = ";"  ) %>% unite(V9, c("b", "a", "c"), sep = "; ") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>%
  subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(gtfall, "/Beatriz_Toledo/nocapMerge/exonscoordinatesallpcstr2filtereddegradationama.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

#on shell do bedtools intersect -a exonscoordinatesallpcstr2filtereddegradationama.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > mergedoverlatirpcstr2filtereddegradationNODES2022.merged.txt

notmergedoverlatirwhip <- read.delim("/Beatriz_Toledo/nocapMerge/mergedoverlatirpcstr2filtereddegradationNODES2022.merged.txt", header=FALSE)
notmergedoverlatirwhipname <- notmergedoverlatirwhip %>%subset(select = c("V10")) %>% dplyr::rename(name=V10)
all.notmerged.simplenamecheck4<-  NODES %>% subset(select = c("name")) %>% distinct()%>% setdiff(whipCEnotmergedname) %>% setdiff(whipAAADnotmergedname) %>% setdiff(whipAAAD1notmergedname) %>%setdiff(notmergedoverlatirwhipname) %>%merge(NODES, by="name") 
head(all.notmerged.simplenamecheck4, n = 20L) 
notmergedoverlatirwhip1<- notmergedoverlatirwhip %>% subset(select = c("V1", "V4","V6", "V8", "V9", "V10")) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1", "V4", "V6", "V8", "V9", "V10"), new = c('Chr','V9', "Strand", "Start", "End", "name"))

EVENTS <- whipCEnotmerged%>% rbind(whipAAADnotmerged) %>%rbind(whipAAAD1notmerged) %>%subset( select= c('Chr','V9', "Strand", "Start", "End", "name")) %>%rbind(notmergedoverlatirwhip1) %>% distinct() %>% mutate (V9 = gsub('"', '', V9))
NODES %>% subset(select= c("name")) %>% distinct() %>%nrow()
EVENTSname <- EVENTS %>% subset(select= c("name")) %>% distinct()
#total 2682, 68 missing

all.ensCE <- ensembl_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name"))
#PB.9743-19, PB.9582-18,PB.9316-9, PB.8772-6,PB.2280-4,  PB.12316-82, PB.4866-14 - all ens
all.ensstart <- ensembl_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(all.ensCE) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End)
#PB.6167-3, PB.5031-60, PB.4800-6, PB.4866-7, PB.2781-9, PB.8482-32 - all ens
all.ensend <- ensembl_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V5"), by.y = c("End"))  %>% anti_join(all.ensCE)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start)
#PB.10632-10, PB.2400-2, PB.10531-13 
all.ensend1<- all.ensend[-grep('PB.4732-34', all.ensend$name), ]

all.gencCE <- gencode_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name"))
all.gencstart <- gencode_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(all.ensCE) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End)
all.gencend <- gencode_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V5"), by.y = c("End"))  %>% anti_join(all.ensCE)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start)

all.RefSeqCE <- ncbiRefSeq2_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4","V5"), by.y = c("Start", "End")) %>%subset(select = c("V1", "V4","V5" ,"V7", "V9", "name"))
all.RefSeqCE1<- all.RefSeqCE[grep('PB.4732-34|PB.1991-12', all.RefSeqCE$name), ]
#PB.4732-34,PB.1991-12
all.RefSeqstart <- ncbiRefSeq2_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V4"), by.y = c("Start"))  %>% anti_join(all.ensCE) %>% subset(select = c("V1", "V4","End" ,"V7", "V9", "name")) %>% dplyr::rename(V5=End)
all.RefSeqstart1<- all.RefSeqstart[grep('PB.4052-20|PB.5326-30|PB.12524-23|PB.3724-21|PB.311-28|PB.1400-48|PB.449-3', all.RefSeqstart$name), ]
PB.4052-20|PB.5326-30|PB.12524-23|PB.3724-21|PB.311-28|PB.1400-48|PB.449-3
all.RefSeqend <- ncbiRefSeq2_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>% merge(all.notmerged.simplenamecheck4, by.x = c("V5"), by.y = c("End"))  %>% anti_join(all.ensCE)%>% subset(select = c("V1", "Start","V5" ,"V7", "V9", "name")) %>% dplyr::rename(V4=Start)
PB.10942-29|PB.8463-5|PB.10400-12
all.RefSeqend1<- all.RefSeqend[grep('PB.10942-29|PB.8463-5|PB.10400-12', all.RefSeqend$name), ]

exonsensembl_SJsqanti_correctedjncmorethan3<- ensembl_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(exonsensembl_SJsqanti_correctedjncmorethan3, "/Beatriz_Toledo/nocapMerge/exonscoordinatesensembl.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
exonsgencode_SJsqanti_correctedjncmorethan3<- gencode_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(exonsgencode_SJsqanti_correctedjncmorethan3, "/Beatriz_Toledo/nocapMerge/exonscoordinatesgenc.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
exonsncbiRefSeq2_SJsqanti_correctedjncmorethan3<- ncbiRefSeq2_SJsqanti_correctedjncmorethan3 %>% subset(V3 == "exon") %>%  mutate(Chr=gsub('^', 'chr', V1))%>% mutate(Zero = 0) %>% subset(V3=="exon", select = c("Chr", "V4", "V5", "V9", "Zero", "V7"))
write.table(exonsncbiRefSeq2_SJsqanti_correctedjncmorethan3, "/Beatriz_Toledo/nocapMerge/exonscoordinatesrefseq.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  

#bedtools intersect -a exonscoordinatesensembl.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > NODES2022.merged.overlatirensemblFilteredincluded.txt
#bedtools intersect -a exonscoordinatesgenc.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > NODES2022.merged.overlatirgencFilteredincluded.txt
#bedtools intersect -a exonscoordinatesrefseq.txt -b NODES2022.merged.eventstocheckIR.txt -wo -s -F 0.99 -filenames > NODES2022.merged.overlatirrefseqFilteredincluded.txt
notmergedoverlatirensembl <- read.delim("/Beatriz_Toledo/nocapMerge/NODES2022.merged.overlatirensemblFilteredincluded.txt", header=FALSE)
notmergedoverlatirgenc <- read.delim("/Beatriz_Toledo/nocapMerge/NODES2022.merged.overlatirgencFilteredincluded.txt", header=FALSE)
notmergedoverlatirrefseq <- read.delim("/Beatriz_Toledo/nocapMerge/NODES2022.merged.overlatirrefseqFilteredincluded.txt", header=FALSE)
genctocheck <- merge(all.notmerged.simplenamecheck4, notmergedoverlatirgenc, by.x = "name", by.y = "V10") 
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
cdscoord <- mergedpcstr2.nofus_NOCAP_sqanti_corrected.gtf.cds %>% subset(V3 == "CDS") %>% separate(V9, into = c("V9", "V10"), sep = ";") %>% mutate (V10 = gsub(' transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10))
cdscoord1 <- cdscoord %>% dplyr::group_by(V10) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V9", "V10","V7", "start", "end")) %>% distinct()
cdscoordens <- Ensembl_corrected.gtf.cds %>% subset(V3 == "CDS") %>% separate(V9, into = c("V9", "V10", "V11"), sep = ";") %>% mutate (V10 = gsub(' transcript_id ', '', V10)) %>% mutate (V10 = gsub('"', '', V10))
cdscoordens1 <- cdscoordens %>% dplyr::group_by(V10) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V9", "V10","V7", "start", "end")) %>% distinct()
cdscoordref <- mm10.ncbiRefSeq2 %>% subset(V3 == "CDS")  %>% separate(V9, into = c("V10", "V11", "V12"), sep = "; ", remove = FALSE) %>% mutate (V11 = gsub('transcript_id ', '', V11)) %>% mutate (V11 = gsub('"', '', V11)) 
cdscoordref1 <- cdscoordref %>% dplyr::group_by(V11) %>% mutate(start = min(V4)) %>% mutate(end = max(V5)) %>% subset(select = c("V10", "V11","V7", "start", "end")) %>% distinct() %>% rename(V9 = V10) %>% rename(V10 = V11)
cdscoordall<- cdscoord1 %>% rbind(cdscoordens1) %>% rbind(cdscoordref1)
EVENTSx<-EVENTS[grep('; tr', EVENTS$V9), ] %>% separate(V9, into = c("gene", "tr", "x"), sep = ";") %>% mutate (tr = gsub(' transcript_id ', '', tr))
EVENTSx2<-EVENTS[-grep('; tr', EVENTS$V9), ] %>% separate(V9, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr))
ANOThave2x2<-ANOThave2 %>% separate(V9, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr))  %>% setnames( old = c("V1", "V7", "V4", "V5"), new = c("Chr","Strand", "Start", "End"))
ANOThave3x2<-ANOThave3 %>% separate(V4, into = c("tr", "gene", "x"), sep = ";") %>% mutate (tr = gsub('transcript_id ', '', tr)) %>% mutate(V1 = gsub("chr", "", V1)) %>% setnames( old = c("V1"), new = c('Chr'))
EVENTSX <- rbind(EVENTSx,EVENTSx2)%>% rbind(ANOThave2x2) %>% rbind(ANOThave3x2)%>% mutate (tr = gsub('"', '', tr))
EVENTSXnodes <- EVENTSX %>% subset(select =c("name"))%>% distinct()
NODESTOBEANALYZED <- intersect(nodesisoexcnodes,EVENTSXnodes)
#2581 NODES FOR NEXT ANALYSES, But have to recheck removing the fusion opnes

mergedpcstr.fusion <- mergedpcstr2.nofus_NOCAP_sqanti_classification%>% subset(structural_category == "fusion",select  = c("isoform")) 
classgeneINCL <- rbind(ensembl_SJsqanti_classification, refseq_SJsqanti_classification) %>%rbind(mergedpcstr2.nofus_NOCAP_sqanti_classification)%>% subset (select =c("isoform", "associated_gene")) %>% merge(mart_exportnoncoding.codingGENE, by.x = "associated_gene", by.y = "Gene.stable.ID", all.x = TRUE)%>%
  merge(EVENTSX, by.x = "isoform", by.y = "tr", all.y = TRUE)  %>%   merge(NODESTOBEANALYZED, by = "name") %>% anti_join(mergedpcstr.fusion)
classgeneINCL<- classgeneINCL%>% mutate(associated_gene=gsub('Sugp2', "ENSMUSG00000036054", associated_gene))%>%mutate(associated_gene=gsub('Brwd1', "ENSMUSG00000022914", associated_gene))%>%mutate(associated_gene=gsub('Rnf214', "ENSMUSG00000042790", associated_gene))%>%mutate(associated_gene=gsub('Gm11266', "ENSMUSG00000087413", associated_gene))%>%
  mutate(associated_gene=gsub('Ubxn8', "ENSMUSG00000052906", associated_gene))%>%mutate(associated_gene=gsub('Ttc28', "ENSMUSG00000033209", associated_gene))%>%mutate(associated_gene=gsub('Poli', "ENSMUSG00000038425", associated_gene))%>%mutate(associated_gene=gsub('Plekha6', "ENSMUSG00000041757", associated_gene))%>%
  mutate(associated_gene=gsub('Med12l', "ENSMUSG00000056476", associated_gene))%>%mutate(associated_gene=gsub('Trip12', "ENSMUSG00000026219", associated_gene))%>%mutate(associated_gene=gsub('Shank1', "ENSMUSG00000038738", associated_gene))%>%mutate(associated_gene=gsub('Mbd6', "ENSMUSG00000025409", associated_gene))%>%
  mutate(associated_gene=gsub('Ksr1', "ENSMUSG00000018334", associated_gene))%>%mutate(associated_gene=gsub('Gm29695', "ENSMUSG00000118219", associated_gene))%>%mutate(associated_gene=gsub('Exoc5', "ENSMUSG00000061244", associated_gene))%>%mutate(associated_gene=gsub('Dlg1', "ENSMUSG00000022770", associated_gene))%>%
  mutate(associated_gene=gsub('Cdh6', "ENSMUSG00000039385", associated_gene))%>%mutate(associated_gene=gsub('Cand2', "ENSMUSG00000030319", associated_gene))%>%mutate(associated_gene=gsub('Unc13a', "ENSMUSG00000034799", associated_gene))%>%mutate(associated_gene=gsub('Tipin', "ENSMUSG00000032397", associated_gene))%>%
  mutate(associated_gene=gsub('Socs7', "ENSMUSG00000038485", associated_gene))%>%mutate(associated_gene=gsub('Slc12a2', "ENSMUSG00000024597", associated_gene))%>%mutate(associated_gene=gsub('Sgsm2', "ENSMUSG00000038351", associated_gene))%>%mutate(associated_gene=gsub('Rapgef2', "ENSMUSG00000062232", associated_gene))%>%
  mutate(associated_gene=gsub('Pip4k2a', "ENSMUSG00000026737", associated_gene))%>%mutate(associated_gene=gsub('Mre11a', "ENSMUSG00000031928", associated_gene))%>%mutate(associated_gene=gsub('Morn1', "ENSMUSG00000029049", associated_gene))%>%mutate(associated_gene=gsub('Lmo3', "ENSMUSG00000030226", associated_gene))%>%
  mutate(associated_gene=gsub('Klhdc8a', "ENSMUSG00000042115", associated_gene))%>%mutate(associated_gene=gsub('Kat2b', "ENSMUSG00000000708", associated_gene))%>%mutate(associated_gene=gsub('Herc4', "ENSMUSG00000020064", associated_gene))%>%mutate(associated_gene=gsub('Fdxr', "ENSMUSG00000018861", associated_gene))%>%
  mutate(associated_gene=gsub('Crim1', "ENSMUSG00000024074", associated_gene))%>%mutate(associated_gene=gsub('Cercam', "ENSMUSG00000039787", associated_gene))%>%mutate(associated_gene=gsub('Ccdc163', "ENSMUSG00000028689", associated_gene))%>%mutate(associated_gene=gsub('Cbwd1', "ENSMUSG00000024878", associated_gene))%>%
  mutate(associated_gene=gsub('Camkmt', "ENSMUSG00000071037", associated_gene))%>%mutate(associated_gene=gsub('Bzrap1', "ENSMUSG00000034156", associated_gene))%>%mutate(associated_gene=gsub('Abcd4', "ENSMUSG00000021240", associated_gene)) %>%mutate(associated_gene=gsub('Ddr1', "ENSMUSG00000003534", associated_gene)) 
classgeneINCLnames <- classgeneINCL %>%subset(select =c("name")) %>% distinct()

NODESTOBEANALYZED2 <- intersect(classgeneINCLnames,classgeneEXCnames)
#2573 NODES TO BE ANALYZED
#Add start and end of cds
classgeneINCL %>% inner_join(NODESTOBEANALYZED2)%>% subset(select =c("associated_gene")) %>% distinct()
#1674 genes
NODESTOBEANALYZED2LEILAPPDP <- rbind(PP.DP.down.2022.merged,PP.DP.up.2022.merged)%>% inner_join(NODESTOBEANALYZED2)
GENESNODESTOBEANALYZED2LEILAPPDP<- classgeneINCL %>% inner_join(NODESTOBEANALYZED2LEILAPPDP)%>% subset(select =c("associated_gene", "name", "isoform")) %>% distinct()
518
NODESTOBEANALYZED2LEILADPN <- rbind(DP.N.down.2022.merged,DP.N.up.2022.merged)%>% inner_join(NODESTOBEANALYZED2)
GENESNODESTOBEANALYZED2LEILADPN<- classgeneINCL %>% inner_join(NODESTOBEANALYZED2LEILADPN)%>% subset(select =c("associated_gene", "name", "isoform")) %>% distinct()
1463
inner_join(NODESTOBEANALYZED2LEILAPPDP,NODESTOBEANALYZED2LEILADPN)
330
inner_join(NODESTOBEANALYZED2LEILAPPDP,NODESTOBEANALYZED2LEILADPN) %>% inner_join(classgeneINCL)%>% subset(select =c("associated_gene")) %>% distinct()
257
