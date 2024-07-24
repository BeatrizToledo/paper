#!/bin/bash

#script for long-read transcripts alignment to the genome

#load module
module load apps/desalt/1.5.5
module load apps/samtools/1.9

#make directory for alignment results
mkdir DESALT

#align high-quality long-reads to the index. Use pre-built mouse RDBG index
deSALT aln /../mouse_index/ /../long_reads_fasta/hq_transcripts.fasta --thread 15 --seed-step 3 --min-chain-score 27 --max-intron-len 400000  -x ccs -O6,24 -M4 -o /../DESALT/highquality_desalt.sam

#sort with samtools
samtools sort -@ 40 -o /../DESALT/highquality_desalt_sorted.sam /../DESALT/highquality_desalt.sam
  
