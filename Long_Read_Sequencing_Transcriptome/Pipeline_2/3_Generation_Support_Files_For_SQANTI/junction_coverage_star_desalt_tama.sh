#!/bin/bash

#script for generating junction coverage files for SQANTI3
#use with the collapsed and filtered file highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf from step 2 

#load modules
module load apps/STAR/2.7.3a

#generate index from desalt+tama gtf file
STAR --runThreadN 65 --runMode genomeGenerate --genomeDir /../SQANTI_INPUT_DESALT_TAMA/ --genomeFastaFiles /../TAMA/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf --sjdbOverhang 99

#calculate junction coverage for each cell type from short-read sequencing data using STAR. Calculate expression for each replicate separately
#will generate multiple files. SJ.out.tab file will be used for SQANTI3
#sample1 (NSC)  - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L222_R1.fastq /../short-reads/L222_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L222_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L355_R1.fastq /../short-reads/L355_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L355_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L393_R1.fastq /../short-reads/L393_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L393_twopassDT --outSAMtype None --twopassMode Basic
#sample2 (NP) - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L223_R2.fastq /../short-reads/L223_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L223_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L356_R2.fastq /../short-reads/L356_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L356_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L394_R2.fastq /../short-reads/L394_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L394_twopassDT --outSAMtype None --twopassMode Basic
#sample3 (N) - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L224_R2.fastq /../short-reads/L224_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L224_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L357_R2.fastq /../short-reads/L357_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L357_twopassDT --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_DESALT_TAMA/  --readFilesIn /../short-reads/L395_R2.fastq /../short-reads/L395_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_DESALT_TAMA/L395_twopassDT --outSAMtype None --twopassMode Basic
