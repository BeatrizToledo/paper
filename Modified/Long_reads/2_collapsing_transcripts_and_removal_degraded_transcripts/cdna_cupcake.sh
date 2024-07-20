#!/bin/bash

# use with the alignment file from gmap
# cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container


#load modules
module load apps/samtools/1.9
module load apps/singularity/3.7.1 

# sort .sam for cDNA Cupcake collapse
samtools sort -o highquality_gmap1_sorted.sam -O sam highquality_gmap1.sam


# cDNA_CUPCAKE collapse
singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif collapse_isoforms_by_sam.py --input hq_transcripts.fasta  -s  highquality_gmap1_sort_cup.sam -o highquality_gmap1_sort_cup --dun-merge-5-shorter --max_fuzzy_junction=5

#filter for degradation
singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif filter_away_subset.py highquality_gmap1_sort_cup.collapsed
