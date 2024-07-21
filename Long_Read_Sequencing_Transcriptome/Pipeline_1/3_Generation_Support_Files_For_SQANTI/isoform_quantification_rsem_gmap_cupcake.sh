#!/bin/bash

#use with the collapsed and filtered file highquality_gmap_sort_cup.collapsed.filtered.gff from the filter_away_subset.py of step 2 

#load modules
module load apps/rsem/ # version?

#prepare the reference
rsem-prepare-reference --num-threads 40 --gtf /../CUPCAKE/highquality_gmap_sort_cup.collapsed.filtered.gff --star Mus_musculus.GRCm38.dna.primary_assembly.fa /../SQANTI_INPUT/GC.rsemSTAR

#calculate isoform expression from short-reads data from the same cell types as for the long-read sequencing,using rsem. Calculate expression for each replicate separately

#sample1 (NSC)  - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L222_R1.fastq /../short-reads/L222_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L222
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L355_R1.fastq /../short-reads/L355_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L355
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L393_R1.fastq /../short-reads/L393_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L393
#sample2 (NP) - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L223_R1.fastq /../short-reads/L223_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L223
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L356_R1.fastq /../short-reads/L356_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L356
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L394_R1.fastq /../short-reads/L394_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L394
#sample3 (N) - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L224_R1.fastq /../short-reads/L224_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L224
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L357_R1.fastq /../short-reads/L357_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L357
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output /../short-reads/L395_R1.fastq /../short-reads/L395_R2.fastq /../SQANTI_INPUT/GC.rsemSTAR /../SQANTI_INPUT/GC.rsemSTAR_L395
