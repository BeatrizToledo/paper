#!/bin/bash

#load module
module load apps/desalt #version?
module load apps/samtools/1.9

#make directory for alignment results
mkdir DESALT

#align high-quality long-reads to the index. Use pre-built mouse RDBG index
deSALT aln ../mouse/ ../long_reads_fasta/hq_transcripts.fasta --thread 15 --seed-step 3 --min-chain-score 27 --max-intron-len 400000  -x ccs -O6,24 -M4 -o ../DESALT/highquality_desalt.sam

#sort with samtools
samtools sort -@ 40 -o /../highquality_desalt.sorted.sam /../highquality_desalt.sam
  
