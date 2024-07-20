#!/bin/bash

# use with the collapsed file from gmap+cupcake highquality_gmap1_sort_cup.collapsed
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 

# get aboundance of the long reads
singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif get_abundance_post_collapse.py highquality_gmap1_sort_cup.collapsed unpolished.cluster_report.csv
