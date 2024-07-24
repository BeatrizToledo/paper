#!/usr/bin/Rscript

#R code for modifying table
#generate an file that associates m…/../ with PB. and with collumns id length is_fl stat pbid
#m64012_190727_053041/105120134/ccs   NA   Y   unique   PB.16.1
#use with filtered file highquality_desalt_sorted_tama.collapsed_nocap.filtered.gtf to obtain the list of transcripts after filtering
# 1       PBRI  transcript        3199736 3202443 .       -       .       gene_id "G1.1"; transcript_id "transcript/132038"; uniq_trans_id "transcript/132038";
#use with highquality_desalt_sorted_tama.collapsed_nocap_read_support.txt generated from tama_read_support_levels.py command
#merge_gene_id   merge_trans_id  gene_read_count trans_read_count        source_line     support_line
#  G1      G1.1    521     3       1       1:m64012_181221_231243/81005164/ccs,m64012_181221_231243/98568367/ccs,m64012_181221_231243/6686664/ccs

library(tibble)
library(plyr)
library(dplyr)
library(tidyr)

read_support <- read.delim("/../TAMA/highquality_desalt_sorted_tama.collapsed_nocap_read_support.txt")

read_support %>% separate(support_line, into = c("b", "support_line"), sep = '";')
  separate_rows(support_line, sep = ",") %>%
  subset(select = c("support_line", "merge_trans_id")) %>%
  add_column(length = "NA", is_fl = "Y", stat = "unique") %>%
  subset(select = c("support_line", "length", "is_fl", "stat", "merge_trans_id")) %>%
  rename(replace = c("support_line" = "id", "merge_trans_id" = "pbid")) %>%  mutate(pbid=gsub('G', '"PB.', pbid))%>%
 write.table("/../TAMA/highquality_desalt_sorted_tama.collapsed_nocap_read_support_modified.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



