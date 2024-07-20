#!/bin/bash

#load the modules
module load apps/STAR/2.7.3a
module load apps/samtools/1.9.

#build index from the merged gtf file
STAR --runMode genomeGenerate --genomeDir /../STAR_align_merged_gtf/ --genomeFastaFiles  /../Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /../mergedpcstr2filtereddegradation.nofus_NOCAP_sqanti_corrected_nofus.gtf   --runThreadN 4

#align short-read data for each replicate
#sample1 (NSC) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L222 --readFilesIn /../fastq/L222_rfp_minus_Track_3438_3447_R1.fastq.gz /../fastq/L222_rfp_minus_Track-3438-3447_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L355 --readFilesIn /../fastq/L355_rfp_minus_Track_3441_3450_R1.fastq.gz /../fastq/L355_rfp_minus_Track-3441-3450_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L393 --readFilesIn /../fastq/L393_rfp_minus_Track_3444_3453_R1.fastq.gz /../fastq/L393_rfp_minus_Track-3444-3453_R2.fastq.gz
#sample2 (NP) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L223 --readFilesIn /../fastq/L223_rfp_plus_Track_3439_3448_R1.fastq.gz /../fastq/L223_rfp_plus_Track-3439-3448_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L356 --readFilesIn /../fastq/L356_rfp_plus_Track_3442_3451_R1.fastq.gz /../fastq/L356_rfp_plus_Track-3442-3451_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome  --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L394 --readFilesIn /../fastq/L394_rfp_plus_Track_3445_3454_R1.fastq.gz /../fastq/L394_rfp_plus_Track-3445-3454_R2.fastq.gz
#sample3 (N) - 3 replicates
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L224 --readFilesIn /../fastq/L224_gfp_plus_Track_3440_3449_R1.fastq.gz /../fastq/L224_gfp_plus_Track-3440-3449_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L357 --readFilesIn /../fastq/L357_gfp_plus_Track_3443_3452_R1.fastq.gz /../fastq/L357_gfp_plus_Track-3443-3452_R2.fastq.gz
STAR --runThreadN 16 --genomeDir /../STAR_align_merged_gtf/genome   --readFilesCommand zcat  --quantMode GeneCounts --outSAMtype BAM Unsorted --outFileNamePrefix L395 --readFilesIn /../fastq/L395_gfp_plus_Track_3446_3455_R1.fastq.gz /../fastq/L395_gfp_plus_Track-3446-3455_R2.fastq.gz

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
