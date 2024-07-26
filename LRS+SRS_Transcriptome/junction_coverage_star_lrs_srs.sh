#!/bin/bash

#script for generating junction coverage files for SQANTI3
#use with the collapsed and filtered file highquality_gmap_sorted_cupcake.collapsed.filtered.gff from step 2 

#load modules
module load apps/STAR/2.7.3a

#generate index from gmap+cupcake gtf file
STAR --runThreadN 65 --runMode genomeGenerate --genomeDir /../MERGED_TRANSCRIPTOME/ --genomeFastaFiles /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /../MERGED_TRANSCRIPTOME/LRS_SRS.filtered.gtf --sjdbOverhang 99

#calculate junction coverage for each cell type from short-read sequencing data using STAR. Calculate expression for each replicate separately

#sample1 (NSC)  - 3 replicates
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L222_R1.fastq  /../short-reads/L222_R2.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L222_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L355_R1.fastq  /../short-reads/L355_R2.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L355_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L393_R1.fastq  /../short-reads/L393_R2.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L393_lrs_srs_sj --outSAMtype None --twopassMode Basic
#sample2 (NP) - 3 replicates
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L223_R2.fastq  /../short-reads/L223_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L223_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L356_R2.fastq  /../short-reads/L356_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L356_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L394_R2.fastq  /../short-reads/L394_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L394_lrs_srs_sj --outSAMtype None --twopassMode Basic
#sample3 (N) - 3 replicates
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L224_R2.fastq  /../short-reads/L224_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L224_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L357_R2.fastq  /../short-reads/L357_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L357_lrs_srs_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../MERGED_TRANSCRIPTOME/ --readFilesIn  /../short-reads/L395_R2.fastq  /../short-reads/L395_R1.fastq --outFileNamePrefix /../MERGED_TRANSCRIPTOME/L395_lrs_srs_sj --outSAMtype None --twopassMode Basic
