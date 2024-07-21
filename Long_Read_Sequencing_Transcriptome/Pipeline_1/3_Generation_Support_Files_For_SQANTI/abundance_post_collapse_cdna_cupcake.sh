#!/bin/bash

#use with the collapsed file from gmap+cupcake highquality_gmap_sorted_cup.collapsed
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 

#get aboundance of the long reads
singularity exec -B /../SQANTI_INPUT/GMAP_CUPCAKE /projects/globalscratch/cdna_cupcake.sif get_abundance_post_collapse.py /../CUPCAKE/highquality_gmap_sorted_cup.collapsed /../long_reads_fasta/unpolished.cluster_report.csv
