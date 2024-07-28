
#!/bin/bash

#load modules
module load apps/hisat2/2.0.4
module load apps/samtools/1.9
module load apps/stringtie/1.3.5


#build index from a reference gtf file using hisat2

# extract splice sites
hisat2_extract_splice_sites.py /../mm10_filtered.gtf > /../MY_INDEX/mm10_filtered.genome.ss

# extract exons
hisat2_extract_exons.py /../mm10_filtered.gtf  > /../MY_INDEX/mm10_filtered.genome.exon

# build index
hisat2-build -p 66 --exon /../MY_INDEX/mm10_filtered.genome.exon --ss /../MY_INDEX/mm10_filtered.genome.ss /../Mus_musculus.GRCm38.dna.primary_ASSEMBLY.fa /../MY_INDEX/mm10.genome_tran


#MAPping SHORT_READS data to the index using Hisat2

# create directory to store MAPping results
mkdir MAP

#MAP each replicate
#sample1 (NSC) - 3 replicates
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L222_R1.fastq -2 /../SHORT_READS/L222_R2.fastq -S /../MAP/L222.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L355_R1.fastq -2 /../SHORT_READS/L355_R2.fastq -S /../MAP/L355.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L393_R1.fastq -2 /../SHORT_READS/L393_R2.fastq -S /../MAP/L393.sam
#sample2 (NP) - 3 replicates
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L223_R1.fastq -2 /../SHORT_READS/L223_R2.fastq -S /../MAP/L223.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L356_R1.fastq -2 /../SHORT_READS/L356_R2.fastq -S /../MAP/L356.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L394_R1.fastq -2 /../SHORT_READS/L394_R2.fastq -S /../MAP/L394.sam
#sample3 (N) - 3 replicates
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L224_R1.fastq -2 /../SHORT_READS/L224_R2.fastq -S /../MAP/L224.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L357_R1.fastq -2 /../SHORT_READS/L357_R2.fastq -S /../MAP/L357.sam
hisat2 -p 36 --dta -x /../MY_INDEX/mm10.genome_tran -1 /../SHORT_READS/L395_R1.fastq -2 /../SHORT_READS/L395_R2.fastq -S /../MAP/L395.sam

#sort the resulting sam files with samtools
#sample1 (NSC) - 3 replicates
samtools sort -@ 40 -o /../MAP/L222.bam /../MAP/L222.sam 
samtools sort -@ 40 -o /../MAP/L355.bam /../MAP/L355.sam 
samtools sort -@ 40 -o /../MAP/L393.bam /../MAP/L393.sam
#sample2 (NP) - 3 replicates
samtools sort -@ 40 -o /../MAP/L223.bam /../MAP/L223.sam 
samtools sort -@ 40 -o /../MAP/L356.bam /../MAP/L356.sam 
samtools sort -@ 40 -o /../MAP/L394.bam /../MAP/L394.sam
#sample3 (N) - 3 replicates
samtools sort -@ 40 -o /../MAP/L224.bam /../MAP/L224.sam 
samtools sort -@ 40 -o /../MAP/L357.bam /../MAP/L357.sam 
samtools sort -@ 40 -o /../MAP/L395.bam /../MAP/L395.sam

# remove SAM files
rm MAP/*.sam


# make ASSEMBLY using Stringtie

# create directory to store ASSEMBLY results
mkdir ASSEMBLY

#assemble each replicate
#sample1 (NSC) - 3 replicates
stringtie /../MAP/L222.bam -l L222 -p 16  -o /../ASSEMBLY/L222.gtf
stringtie /../MAP/L355.bam -l L355 -p 16  -o /../ASSEMBLY/L355.gtf
stringtie /../MAP/L393.bam -l L393 -p 16  -o /../ASSEMBLY/L393.gtf
#sample2 (NP) - 3 replicates
stringtie /../MAP/L223.bam -l L223 -p 16  -o /../ASSEMBLY/L223.gtf
stringtie /../MAP/L356.bam -l L356 -p 16  -o /../ASSEMBLY/L356.gtf
stringtie /../MAP/L394.bam -l L394 -p 16  -o /../ASSEMBLY/L394.gtf
#sample3 (N) - 3 replicates
stringtie /../MAP/L224.bam -l L224 -p 16  -o /../ASSEMBLY/L224.gtf
stringtie /../MAP/L357.bam -l L357 -p 16  -o /../ASSEMBLY/L357.gtf
stringtie /../MAP/L395.bam -l L395 -p 16  -o /../ASSEMBLY/L395.gtf

#assemble the replicates gtfs (list them in the mergelist) into a uique ASSEMBLY using a minimim coverage of 10 reads (c10); do not include the reference gtf in the final ASSEMBLY 
stringtie --merge -c10 -p 16  -o /../ASSEMBLY/SRS.gtf /../mergelist.txt
