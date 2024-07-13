#!/bin/bash

# MODIFY ID OF PACBIO High Quality transcripts FASTA FILE
# Fasta and csv files given by the facility had the isoform ID unmatched
# High_QualityIsoforms_1878_hq_transcripts.fasta
# Fasta file transcript ID is UnnamedSample_HQ_transcript/0
# unpolished.cluster_report.csv the transcript ID is transcript/0,m54345U_200310_133837/123864620/ccs,FL
# this difference is constant, therefore is just to cut the "UnnamedSample_HQ_" in the fasta file

cat High_QualityIsoforms_1878_hq_transcripts.fasta | sed 's/UnnamedSample_HQ_//g' > hq_transcripts_renamed.fasta

#desalt

./deSALT aln ../mouse/ ../long_reads_fasta/hq_transcripts.fasta --thread 15 --seed-step 3 --min-chain-score 27 --max-intron-len 400000  -x ccs -O6,24 -M4 -o hq.sam
