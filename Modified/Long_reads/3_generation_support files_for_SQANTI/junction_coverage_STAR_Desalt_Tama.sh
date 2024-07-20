#!/bin/bash

#script to generate junction coverage for SQANTI3
#use with output of desalt+tama filterdegradationtama.modified.gtf

#load modules
module load apps/STAR/

#generate index from desalt+tama gtf file
STAR --runThreadN 65 --runMode genomeGenerate --genomeDir /Beatriz_Toledo/ --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile filterdegradationtama.modified.gtf --sjdbOverhang 99


#calculate junction coverage for each cell type from short-read sequencing data using STAR. Calculate expression for each replicate separately

#sample1 (NSC)  - 3 replicates
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L222_R1.fastq /Beatriz_Toledo/sj/L222_R2.fastq --outFileNamePrefix L222_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L355_R1.fastq /Beatriz_Toledo/sj/L355_R2.fastq --outFileNamePrefix L355_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L393_R1.fastq /Beatriz_Toledo/sj/L393_R2.fastq --outFileNamePrefix L393_twopassGC --outSAMtype None --twopassMode Basic
#sample2 (NP) - 3 replicates
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L223_R2.fastq /Beatriz_Toledo/sj/L223_R1.fastq --outFileNamePrefix L223_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L356_R2.fastq /Beatriz_Toledo/sj/L356_R1.fastq --outFileNamePrefix L356_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L394_R2.fastq /Beatriz_Toledo/sj/L394_R1.fastq --outFileNamePrefix L394_twopassGC --outSAMtype None --twopassMode Basic
#sample3 (N) - 3 replicates
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L224_R2.fastq /Beatriz_Toledo/sj/L224_R1.fastq --outFileNamePrefix L224_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L357_R2.fastq /Beatriz_Toledo/sj/L357_R1.fastq --outFileNamePrefix L357_twopassGC --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /Beatriz_Toledo/ --readFilesIn /Beatriz_Toledo/sj/L395_R2.fastq /Beatriz_Toledo/sj/L395_R1.fastq --outFileNamePrefix L395_twopassGC --outSAMtype None --twopassMode Basic
