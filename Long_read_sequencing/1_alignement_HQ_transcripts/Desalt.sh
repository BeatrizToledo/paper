#!/bin/bash

./deSALT aln ../mouse/ ../long_reads_fasta/hq_transcripts.fasta --thread 15 --seed-step 3 --min-chain-score 27 --max-intron-len 400000  -x ccs -O6,24 -M4 -o hq.sam
