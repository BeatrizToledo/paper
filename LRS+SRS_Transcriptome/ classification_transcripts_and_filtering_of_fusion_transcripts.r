#!/usr/bin/env Rscript

#R code for removing fusion trancripts
#should be in folder /../MERGED_TRANSCRIPTOME/

library(dplyr) 

#classification ensembl
LRS_SRS <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_classification.txt", header=FALSE,  quote="'")
LRS_SRS_nofusion_classification <- LRS_SRS %>% subset(structural_category != "fusion", select = c("isoform", "structural_category")) %>% rename_with(.cols = c(2,3), ~glue::glue("ens{.}"))

#classification gencode
LRS_SRS_gencode_class <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_gencode_classification.txt", header=FALSE,  quote="'")
LRS_SRS_gencode <- LRS_SRS_gencode_class %>% subset(structural_category != "fusion",select  = c("isoform", "structural_category")) %>% rename_with(.cols = c(2,3), ~glue::glue("gen{.}"))

#classification refseq
LRS_SRS_refreq_class <- read.delim("/../MERGED_TRANSCRIPTOME/LRS_SRS_gencode_classification.txt", header=FALSE,  quote="'")
LRS_SRS_refseq <- LRS_SRS_refreq_class %>% subset(structural_category != "fusion",select  = c("isoform", "structural_category")) %>% rename_with(.cols = c(2,3), ~glue::glue("ref{.}"))

#table with classification of ensembl, refseq and gencode
classmerge <- LRS_SRS_ensembl %>% merge(LRS_SRS_gencode, by = "isoform", all.x = TRUE) %>% merge(LRS_SRS_refseq, by = "isoform", all.x = TRUE)

#
classmergefsm <- classmerge %>% subset(ensstructural_category == "full-splice_match" | genstructural_category == "full-splice_match" |  refstructural_category == "full-splice_match") 
classmergepsm <- classmerge %>% setdiff(classmergefsm) %>%subset(ensstructural_category == "incomplete-splice_match" | genstructural_category == "incomplete-splice_match" | refstructural_category == "incomplete-splice_match") 
classmergenic <- classmerge %>% setdiff(classmergefsm) %>% setdiff(classmergepsm) %>%subset(ensstructural_category == "novel_in_catalog" | genstructural_category == "novel_in_catalog"  | refstructural_category == "novel_in_catalog")
classmergennc <- classmerge %>% setdiff(classmergefsm) %>%  setdiff(classmergepsm) %>% setdiff(classmergenic) %>%subset(ensstructural_category == "novel_not_in_catalog" | genstructural_category == "novel_not_in_catalog" | refstructural_category == "novel_not_in_catalog")
classmergeothers <- classmerge %>% setdiff(classmergefsm) %>%  setdiff(classmergepsm) %>% setdiff(classmergenic) %>% setdiff(classmergennc) 
classmergefsm1 <- classmergefsm %>% subset(select = c("isoform"))%>% mutate(structural_category = "full-splice match")
classmergenic1 <- classmergenic %>% subset(select = c("isoform")) %>% mutate(structural_category = "novel in catalog")
classmergennc1 <- classmergennc %>% subset(select = c("isoform")) %>% mutate(structural_category = "novel not in catalog")
classmergepsm1 <- classmergepsm %>% subset(select = c("isoform")) %>% mutate(structural_category = "partial splice match")
classmergeothers1 <- classmergeothers  %>% subset(select = c("isoform")) %>% mutate(structural_category = "others")

#final table with classification
LRS_SRS_nofusion_classification <- rbind(classmergefsm1,classmergenic1) %>% rbind(classmergennc1)%>% rbind(classmergeism1) %>% rbind(classmergeothers1)
write.table(LRS_SRS_nofusion_classification,"/../MERGED_TRANSCRIPTOME/LRS_SRS_nofusion_classification.txt", sep="\t", quote= FALSE , row.name=FALSE, col.names = FALSE)
