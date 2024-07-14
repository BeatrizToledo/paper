#!/bin/bash

# I need to generate FL count for tama transcripts, to use it as input of SQANTI
# in R 

Rscript tama_FL.r

# In the command line in the terminal

cat filterdegradationtama_nocap_trans_readsupport_modified.txt | sed 's/tama://g' > filterdegradationtama_nocap_trans_readsupport_modified2.txt
cat filterdegradationtama_nocap_trans_readsupport_modified2.txt | sed 's/G/PB./g' > filterdegradationtama_nocap_trans_readsupport_modified3.txt
python fa2fq.py filterdegradationtama.modified.fa
python demux_isoseq_with_genome.py --mapped_fafq  filterdegradationtama.modified.fastq --read_stat filterdegradationtama_nocap_trans_readsupport_modified3.txt --classify_csv flnc.report.csv -o filterdegradationtama.mapped_fl_count.txt
