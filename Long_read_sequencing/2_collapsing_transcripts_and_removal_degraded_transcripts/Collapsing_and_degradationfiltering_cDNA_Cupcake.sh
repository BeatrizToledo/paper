#!/bin/bash

# SORT SAM for cDNA Cupcake collapse
# Optional -l 0 - to not compress -@ 10 
# Different way of sorting: sort -k 3,3 -k 4,4n highquality_gmap1.sam > highquality_gmap1_sort_cup.sam

samtools sort -o highquality_gmap1_sorted.sam -O sam highquality_gmap1.sam


# cDNA_CUPCAKE collapse

module load apps/singularity/3.7.1 
singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif collapse_isoforms_by_sam.py --input hq_transcripts_renamed.fasta  -s  highquality_gmap1_sort_cup.sam -o highquality_gmap1_sort_cup --dun-merge-5-shorter --max_fuzzy_junction=5

#filter for degradation

singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif filter_away_subset.py highquality_gmap1_sort_cup.collapsed
