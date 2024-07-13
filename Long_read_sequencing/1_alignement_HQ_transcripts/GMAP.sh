#!/bin/bash

# MODIFY ID OF PACBIO High Quality transcripts FASTA FILE
# Fasta and csv files given by the facility had the isoform ID unmatched
# High_QualityIsoforms_1878_hq_transcripts.fasta
# Fasta file transcript ID is UnnamedSample_HQ_transcript/0
# unpolished.cluster_report.csv the transcript ID is transcript/0,m54345U_200310_133837/123864620/ccs,FL
# this difference is constant, therefore is just to cut the "UnnamedSample_HQ_" in the fasta file

cat High_QualityIsoforms_1878_hq_transcripts.fasta | sed 's/UnnamedSample_HQ_//g' > hq_transcripts_renamed.fasta

#Kept the same size, same number of lines, and the head is  >transcript/0

# now ALIGN USING GMAP
# -n 1 - just 1 output
# --split-output=STRING          Basename for multiple-file output, separately for nomapping, 
#                                   uniq, mult, (and chimera, if --chimera-margin is selected)

gmap_build -d GENOME_Assembly -s none -k 15 -D /Long-read_sequencing/GMAP/ Mus_musculus.GRCm38.dna.primary_assembly.fa
gmap -D /GMAP/GMAP/ -d GENOME_Assembly hq_transcripts_renamed.fasta --no-chimeras --min-trimmed-coverage 0.85 --min-identity 0.9 -B 5 -t 24 -f samse --cross-species -z sense_force -n 1 -K 400000 > highquality_gmap1.sam
