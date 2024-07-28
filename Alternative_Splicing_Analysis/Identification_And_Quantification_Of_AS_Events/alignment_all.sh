#!/bin/bash

#code to generate aligment file for Whippet, to allow for the discovery of cryptic events

#load the modules
module load apps/STAR/2.7.3a
module load apps/samtools/1.9

#build index from the merged gtf file
STAR --runMode genomeGenerate --genomeDir /../STAR_ALIGN_MERGED_GTF/ --genomeFastaFiles  /../Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /../mergedpcstr2filtereddegradation.nofus_NOCAP_sqanti_corrected_nofus.gtf   --runThreadN 4

#align short-read paired-end sequencing data (_R1.fastq.gz and _R2.fastq.gz) for each replicate 
#sample1 (NSC) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L222 --readFilesIn /../SHORT_READS/L222_R1.fastq.gz /../SHORT_READS/L222_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L355 --readFilesIn /../SHORT_READS/L355_R1.fastq.gz /../SHORT_READS/L355_R2.fastq.gz 
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L393 --readFilesIn /../SHORT_READS/L393_R1.fastq.gz /../SHORT_READS/L393_R2.fastq.gz
#sample2 (NP) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L223 --readFilesIn /../SHORT_READS/L223_R1.fastq.gz /../SHORT_READS/L223_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L356 --readFilesIn /../SHORT_READS/L356_R1.fastq.gz /../SHORT_READS/L356_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L394 --readFilesIn /../SHORT_READS/L394_R1.fastq.gz /../SHORT_READS/L394_R2.fastq.gz
#sample3 (N) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L224 --readFilesIn /../SHORT_READS/L224_R1.fastq.gz /../SHORT_READS/L224_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L357 --readFilesIn /../SHORT_READS/L357_R1.fastq.gz /../SHORT_READS/L357_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_ALIGN_MERGED_GTF/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L395 --readFilesIn /../SHORT_READS/L395_R1.fastq.gz /../SHORT_READS/L395_R2.fastq.gz

#merge single .bam files into a unique one with samtools merge
samtools merge All.merged.bam \
    L222Aligned.out.bam \
    L355Aligned.out.bam \
    L393Aligned.out.bam \
    L223Aligned.out.bam \
    L356Aligned.out.bam \
    L394Aligned.out.bam \
    L224Aligned.out.bam \
    L357Aligned.out.bam \
    L395Aligned.out.bam


#sort the merged .bam file with samtools sort
samtools sort -o All.merged.bam All.merged.sorted

#remove duplicates with samtools rmdup
samtools rmdup -S All.merged.sorted All.merged.sorted.rmdup.bam

#create an index for the resulting bam file
samtools index All.merged.sorted.rmdup.bam
