#!/bin/bash

#use with the alignment file from gmap highquality_gmap.sam (from step 1)
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 
module load apps/cupcake/12.0.0

#make directory for output
mkdir CUPCAKE

#cdna_cupcake collapse
singularity exec -B /../CUPCAKE /projects/globalscratch/cdna_cupcake.sif collapse_isoforms_by_sam.py --input /../long_reads_fasta/hq_transcripts.fasta  -s  /../GMAP/highquality_gmap_sorted.sam -o  /../CUPCAKE/highquality_gmap_sorted_cup --dun-merge-5-shorter --max_fuzzy_junction=5

#filter for degradation
singularity exec -B /../CUPCAKE /projects/globalscratch/cdna_cupcake.sif filter_away_subset.py /../CUPCAKE/highquality_gmap_sorted_cup.collapsed

#make directory for resulst of step 3
mkdir SQANTI_INPUT_GMAP_CUPCAKE
