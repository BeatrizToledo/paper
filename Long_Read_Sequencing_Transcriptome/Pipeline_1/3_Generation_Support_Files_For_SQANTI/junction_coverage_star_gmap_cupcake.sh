#!/bin/bash

#script for generating junction coverage files for SQANTI3
#use with the collapsed and filtered file highquality_gmap_sorted_cupcake.collapsed.filtered.gff from step 2 

#load modules
module load apps/STAR/2.7.3a

#generate index from gmap+cupcake gtf file
STAR --runThreadN 65 --runMode genomeGenerate --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --genomeFastaFiles /../CUPCAKE/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.gff --sjdbOverhang 99

#calculate junction coverage for each cell type from short-read sequencing data using STAR. Calculate expression for each replicate separately

#sample1 (NSC)  - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L222_R1.fastq  /../short-reads/L222_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L222_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L355_R1.fastq  /../short-reads/L355_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L355_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L393_R1.fastq  /../short-reads/L393_R2.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L393_gc_sj --outSAMtype None --twopassMode Basic
#sample2 (NP) - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L223_R2.fastq  /../short-reads/L223_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L223_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L356_R2.fastq  /../short-reads/L356_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L356_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L394_R2.fastq  /../short-reads/L394_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L394_gc_sj --outSAMtype None --twopassMode Basic
#sample3 (N) - 3 replicates
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L224_R2.fastq  /../short-reads/L224_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L224_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L357_R2.fastq  /../short-reads/L357_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L357_gc_sj --outSAMtype None --twopassMode Basic
STAR --runThreadN 65 --genomeDir /../SQANTI_INPUT_GMAP_CUPCAKE/ --readFilesIn  /../short-reads/L395_R2.fastq  /../short-reads/L395_R1.fastq --outFileNamePrefix /../SQANTI_INPUT_GMAP_CUPCAKE/L395_gc_sj --outSAMtype None --twopassMode Basic
