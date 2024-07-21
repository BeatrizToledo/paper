
#!/bin/bash

#load modules
module load apps/hisat2/2.0.4
module load apps/samtools/1.9
module load apps/stringtie/1.3.5


#build index from a reference gtf file using hisat2

# extract splice sites
hisat2_extract_splice_sites.py /../mm10_filtered.gtf > /../my_index/mm10_filtered.genome.ss

# extract exons
hisat2_extract_exons.py /../mm10_filtered.gtf  > /../my_index/mm10_filtered.genome.exon

# build index
hisat2-build -p 66 --exon /../my_index/mm10_filtered.genome.exon --ss /../my_index/mm10_filtered.genome.ss /../Mus_musculus.GRCm38.dna.primary_assembly.fa /../my_index/mm10.genome_tran


#mapping short-reads data to the index using Hisat2

# create directory to store mapping results
mkdir map

#map each replicate
#sample1 (NSC) - 3 replicates
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L222_R1.fastq -2 /../short-reads/L222_R2.fastq -S /../map/L222.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L355_R1.fastq -2 /../short-reads/L355_R2.fastq -S /../map/L355.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L393_R1.fastq -2 /../short-reads/L393_R2.fastq -S /../map/L393.sam
#sample2 (NP) - 3 replicates
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L223_R1.fastq -2 /../short-reads/L223_R2.fastq -S /../map/L223.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L356_R1.fastq -2 /../short-reads/L356_R2.fastq -S /../map/L356.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L394_R1.fastq -2 /../short-reads/L394_R2.fastq -S /../map/L394.sam
#sample3 (N) - 3 replicates
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L224_R1.fastq -2 /../short-reads/L224_R2.fastq -S /../map/L224.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L357_R1.fastq -2 /../short-reads/L357_R2.fastq -S /../map/L357.sam
hisat2 -p 36 --dta -x /../my_index/mm10.genome_tran -1 /../short-reads/L395_R1.fastq -2 /../short-reads/L395_R2.fastq -S /../map/L395.sam

#sort the resulting sam files with samtools
#sample1 (NSC) - 3 replicates
samtools sort -@ 40 -o /../map/L222.bam /../map/L222.sam 
samtools sort -@ 40 -o /../map/L355.bam /../map/L355.sam 
samtools sort -@ 40 -o /../map/L393.bam /../map/L393.sam
#sample2 (NP) - 3 replicates
samtools sort -@ 40 -o /../map/L223.bam /../map/L223.sam 
samtools sort -@ 40 -o /../map/L356.bam /../map/L356.sam 
samtools sort -@ 40 -o /../map/L394.bam /../map/L394.sam
#sample3 (N) - 3 replicates
samtools sort -@ 40 -o /../map/L224.bam /../map/L224.sam 
samtools sort -@ 40 -o /../map/L357.bam /../map/L357.sam 
samtools sort -@ 40 -o /../map/L395.bam /../map/L395.sam

# remove SAM files
rm map/*.sam


# make assembly using Stringtie

# create directory to store assembly results
mkdir assembly

#assemble each replicate
#sample1 (NSC) - 3 replicates
stringtie /../map/L222.bam -l L222 -p 16  -o /../assembly/L222.gtf
stringtie /../map/L355.bam -l L355 -p 16  -o /../assembly/L355.gtf
stringtie /../map/L393.bam -l L393 -p 16  -o /../assembly/L393.gtf
#sample2 (NP) - 3 replicates
stringtie /../map/L223.bam -l L223 -p 16  -o /../assembly/L223.gtf
stringtie /../map/L356.bam -l L356 -p 16  -o /../assembly/L356.gtf
stringtie /../map/L394.bam -l L394 -p 16  -o /../assembly/L394.gtf
#sample3 (N) - 3 replicates
stringtie /../map/L224.bam -l L224 -p 16  -o /../assembly/L224.gtf
stringtie /../map/L357.bam -l L357 -p 16  -o /../assembly/L357.gtf
stringtie /../map/L395.bam -l L395 -p 16  -o /../assembly/L395.gtf

#assemble the replicates gtfs (list them in the mergelist) into a uique assembly using a minimim coverage of 10 reads (c10); do not include the reference gtf in the final assembly 
stringtie --merge -c10 -p 16  -o /../assembly/stringtie.mm10.noG10cov.merged.gtf /../mergelist.txt
