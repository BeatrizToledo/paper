#!/bin/bash

#use with the collapsed and filtered file from gmap+cupcake highquality_gmap1_sort_cup.collapsed
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 
module load apps/rsem/

#prepare the reference
rsem-prepare-reference --num-threads 40 --gtf highquality_gmap1_sort_cup.collapsed.filtered.gff --star Mus_musculus.GRCm38.dna.primary_assembly.fa GC.rsemSTAR

#calculate isoform expression from short-reads data from the same cell types as for the long-read sequencing,using rsem. Calculate expression for each replicate separately

#sample1 (NSC)  - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L222_R1.fastq L222_R2.fastq GC.rsemSTAR GC.rsemSTAR_L222
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L355_R1.fastq L355_R2.fastq GC.rsemSTAR GC.rsemSTAR_L355
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L393_R1.fastq L393_R2.fastq GC.rsemSTAR GC.rsemSTAR_L393
#sample2 (NP) - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L223_R1.fastq L223_R2.fastq GC.rsemSTAR GC.rsemSTAR_L223
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L356_R1.fastq L356_R2.fastq GC.rsemSTAR GC.rsemSTAR_L356
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L394_R1.fastq L394_R2.fastq GC.rsemSTAR GC.rsemSTAR_L394
#sample3 (N) - 3 replicates
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L224_R1.fastq L224_R2.fastq GC.rsemSTAR GC.rsemSTAR_L224
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L357_R1.fastq L357_R2.fastq GC.rsemSTAR GC.rsemSTAR_L357
rsem-calculate-expression --paired-end --strand-specific --star --num-threads 45 --no-bam-output L395_R1.fastq L395_R2.fastq GC.rsemSTAR GC.rsemSTAR_L395
