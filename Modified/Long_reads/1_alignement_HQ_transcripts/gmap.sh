#!/bin/bash

#load modules
module load apps/gmap

#build index assembly from reference fasta file
gmap_build -d GENOME_Assembly -s none -k 15 -D /Long-read_sequencing/GMAP/ Mus_musculus.GRCm38.dna.primary_assembly.fa

##align high-quality long-reads to the index
gmap -D /GMAP/GMAP/ -d GENOME_Assembly hq_transcripts.fasta --no-chimeras --min-trimmed-coverage 0.85 --min-identity 0.9 -B 5 -t 24 -f samse --cross-species -z sense_force -n 1 -K 400000 > highquality_gmap1.sam
