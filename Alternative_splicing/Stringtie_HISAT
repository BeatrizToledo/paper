#!/bin/bash

module load apps/stringtie/1.3.5

#first see with only the ensembl gtf

stringtie /Beatriz_Toledo/Stringtie/L222.bam -l L222 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L222.gtf
stringtie /Beatriz_Toledo/Stringtie/L223.bam -l L223 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L223.gtf
stringtie /Beatriz_Toledo/Stringtie/L224.bam -l L224 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L224.gtf
stringtie /Beatriz_Toledo/Stringtie/L355.bam -l L355 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L355.gtf
stringtie /Beatriz_Toledo/Stringtie/L356.bam -l L356 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L356.gtf
stringtie /Beatriz_Toledo/Stringtie/L357.bam -l L357 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L357.gtf
stringtie /Beatriz_Toledo/Stringtie/L393.bam -l L393 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L393.gtf
stringtie /Beatriz_Toledo/Stringtie/L394.bam -l L394 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L394.gtf
stringtie /Beatriz_Toledo/Stringtie/L395.bam -l L395 -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/L395.gtf

#merge only mm10
stringtie --merge -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.mm10.merged.gtf /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/mergelist.txt
stringtie --merge -c10 -p 16  -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.mm10.noG10cov.merged.gtf /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/mergelist.txt

#merge also with PacBio
stringtie --merge -p 16 -G /Beatriz_Toledo/Stringtie/mm10.filtered.gtf -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.merged.gtf /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/mergelist2.txt
stringtie --merge -c10 -p 16  -o /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/mergelist2.txt
 
#re-qunatify
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L222/L222.gtf /Beatriz_Toledo/Stringtie/L222.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L223/L223.gtf /Beatriz_Toledo/Stringtie/L223.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L224/L224.gtf /Beatriz_Toledo/Stringtie/L224.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L355/L355.gtf /Beatriz_Toledo/Stringtie/L355.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L356/L356.gtf /Beatriz_Toledo/Stringtie/L356.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L357/L357.gtf /Beatriz_Toledo/Stringtie/L357.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L393/L393.gtf /Beatriz_Toledo/Stringtie/L393.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L394/L394.gtf /Beatriz_Toledo/Stringtie/L394.bam
stringtie -e -B -p 16 -G /Beatriz_Toledo/Stringtie/assembly/mm10_ensembl_filtered/stringtie.all.noG10cov.merged.gtf -o /Beatriz_Toledo/Stringtie/ballgown/mm10_ensembl_filtered/L395/L395.gtf /Beatriz_Toledo/Stringtie/L395.bam

module load apps/hisat2/2.0.4
module load apps/samtools/1.9

#Merged stringtie.mm10.noG10cov.merged.strand.fixed.bed + Pacbio.bed = PCST.bed

#build index
hisat2_extract_splice_sites.py /Beatriz_Toledo/Whippet.jl/PCST.gtf > /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome.ss

hisat2_extract_exons.py /Beatriz_Toledo/Whippet.jl/PCST.gtf > /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome.exon

hisat2-build -p 66 --exon /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome.exon --ss /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome.ss /Beatriz_Toledo/Whippet.jl/Mus_musculus.GRCm38.dna.primary_assembly.fa /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran

#mapping

hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L222_R1.fastq -2 /Beatriz_Toledo/Whippet.jl/L222_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L222.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L223_R1.fastq -2 /Beatriz_Toledo/Whippet.jl/L223_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L223.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L224_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L224_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L224.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L355_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L355_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L355.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L356_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L356_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L356.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L357_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L357_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L357.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L393_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L393_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L393.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L394_R1.fastq  -2 /Beatriz_Toledo/Whippet.jl/L394_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L394.sam
hisat2 -p 36 --dta -x /Beatriz_Toledo/Whippet.jl/my_index/mm10.genome_tran -1 /Beatriz_Toledo/Whippet.jl/L395_R1.fastq -2 /Beatriz_Toledo/Whippet.jl/L395_R2.fastq -S /Beatriz_Toledo/Whippet.jl/L395.sam

#sort
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L222.bam /Beatriz_Toledo/Whippet.jl/L222.sam 
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L223.bam /Beatriz_Toledo/Whippet.jl/L223.sam 
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L224.bam /Beatriz_Toledo/Whippet.jl/L224.sam 
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L355.bam /Beatriz_Toledo/Whippet.jl/L355.sam 
samtools sort -@ 60 -o /Beatriz_Toledo/Whippet.jl/L356.bam /Beatriz_Toledo/Whippet.jl/L356.sam 
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L357.bam /Beatriz_Toledo/Whippet.jl/L357.sam 
samtools sort -@ 60 -o /Beatriz_Toledo/Whippet.jl/L393.bam /Beatriz_Toledo/Whippet.jl/L393.sam
samtools sort -@ 40 -o /Beatriz_Toledo/Whippet.jl/L394.bam /Beatriz_Toledo/Whippet.jl/L394.sam
samtools sort -@ 60 -o /Beatriz_Toledo/Whippet.jl/L395.bam /Beatriz_Toledo/Whippet.jl/L395.sam

