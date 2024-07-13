#!/bin/bash

# I need to generate FL count for tama transcripts, to use it as input of SQANTI
# in R 

R

filterdegradationtama_trans_report <- read.delim("/Beatriz_Toledo/tama-master/filterdegradationtama_trans_report.txt")
tama_collapse_cluster_nocap_trans_readsupport_read_support <- read.delim("/Beatriz_Toledo/tama-master/tama_collapse_cluster_nocap_trans_readsupport_read_support.txt")

look_up_table <- filterdegradationtama_trans_report %>% subset(select = c("transcript_id", "all_source_trans")) %>% mutate(all_source_trans=gsub('filtertama_', "", all_source_trans)) %>% dplyr::rename(G = all_source_trans) %>% separate_rows("G", sep = ",")
tab <- tama_collapse_cluster_nocap_trans_readsupport_read_support %>% subset(select =c("merge_trans_id", "trans_read_count")) 
filterdegradationtama_count<-merge(look_up_table, tab,by.x = "G", by.y="merge_trans_id") %>% subset(select = c("transcript_id", "trans_read_count"))%>% group_by(transcript_id) %>%   summarise(trans_read_count = sum(trans_read_count))

#I need to generate an file that associates m…/../ with PB. I think I need to generate it from the Trans_read. But I will also generate from the .bed
# id length is_fl stat pbid
# m64012_190727_053041/105120134/ccs NA Y unique PB.16.1
#This didn’t work so good because it considers the G.1.1 like a gene and the transcript/… like transcript ID. So I will try to modify this in R. python 

tama_collapse_cluster_nocap_trans_readsupport_read_support %>%
  separate_rows(support_line, sep = ",") %>%
  subset(select = c(support_line, merge_trans_id)) %>%
  add_column(length = "NA", is_fl = "Y", stat = "unique") %>%
  subset(select = c(support_line, length, is_fl, stat, merge_trans_id)) %>% merge(look_up_table, by.x = "merge_trans_id", by.y = "G")%>%
  dplyr::rename("id" = "support_line" ,"pbid"= "merge_trans_id") %>% distinct() %>%
  write.table("/Beatriz_Toledo/tama-master/filterdegradationtama_nocap_trans_readsupport_modified.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

q()

# In the command line in the terminal

cat filterdegradationtama_nocap_trans_readsupport_modified.txt | sed 's/tama://g' > filterdegradationtama_nocap_trans_readsupport_modified2.txt
cat filterdegradationtama_nocap_trans_readsupport_modified2.txt | sed 's/G/PB./g' > filterdegradationtama_nocap_trans_readsupport_modified3.txt
python fa2fq.py filterdegradationtama.modified.fa
python demux_isoseq_with_genome.py --mapped_fafq  filterdegradationtama.modified.fastq --read_stat filterdegradationtama_nocap_trans_readsupport_modified3.txt --classify_csv flnc.report.csv -o filterdegradationtama.mapped_fl_count.txt
