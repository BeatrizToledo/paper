#!/bin/bash

#script for long-read transcripts alignment to the genome
#load modules
module load apps/gmap/2019-06-10
module load apps/samtools/1.9

#make directory for alignment results
mkdir GMAP

#build index assembly from reference fasta file
gmap_build -d GENOME_Assembly -s none -k 15 -D /../GMAP/ /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa

#align high-quality long-reads to the index
gmap -D /../GMAP/ -d GENOME_Assembly /../LONG_READS_ISOSEQ3/hq_transcripts.fasta --no-chimeras --min-trimmed-coverage 0.85 --min-identity 0.9 -B 5 -t 24 -f samse --cross-species -z sense_force -n 1 -K 400000 > /../GMAP/highquality_gmap.sam

#sort .sam file
sort -k 3,3 -k 4,4n /../GMAP/highquality_gmap.sam > /../GMAP/highquality_gmap_sorted.sam
