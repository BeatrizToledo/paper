#!/bin/bash

#script for collapsing redundant transcripts and removing potentially degraded transcripts
#use with the alignment file from gmap highquality_gmap_sorted.sam (from step 1)
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 

#make directory for output
mkdir CUPCAKE

#collapsing of redundant transcripts. collapse_isoforms_by_sam.py command generates multiple files, including .gff file.
singularity exec -B /../CUPCAKE,/../LONG_READS_FASTA,/../GMAP /projects/globalscratch/cdna_cupcake.sif collapse_isoforms_by_sam.py --input /../LONG_READS_FASTA/hq_transcripts.fasta  -s  /../GMAP/highquality_gmap_sorted.sam -o  /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed --dun-merge-5-shorter --max_fuzzy_junction=5

#filter out degraded transcripts. filter_away_subset.py command generates multiple files, including .gff and .abundance.txt 
singularity exec -B /../CUPCAKE /projects/globalscratch/cdna_cupcake.sif filter_away_subset.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed

#make directory for resulst of step 3 and 4
mkdir SQANTI_INPUT_GMAP_CUPCAKE
mkdir SQANTI_GMAP_CUPCAKE
