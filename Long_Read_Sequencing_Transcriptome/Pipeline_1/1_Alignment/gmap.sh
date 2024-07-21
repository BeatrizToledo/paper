#!/bin/bash

#load modules
module load apps/gmap #version?
module load apps/samtools/1.9

#make directory for alignment results
mkdir GMAP

#build index assembly from reference fasta file
gmap_build -d GENOME_Assembly -s none -k 15 -D /../GMAP/ Mus_musculus.GRCm38.dna.primary_assembly.fa

#align high-quality long-reads to the index
gmap -D /../GMAP/ -d GENOME_Assembly /../long_reads_fasta/hq_transcripts.fasta --no-chimeras --min-trimmed-coverage 0.85 --min-identity 0.9 -B 5 -t 24 -f samse --cross-species -z sense_force -n 1 -K 400000 > highquality_gmap.sam

#sort .sam for cDNA Cupcake collapse
samtools sort -o /../GMAP/highquality_gmap_sorted.sam -O sam /../GMAP/highquality_gmap.sam
