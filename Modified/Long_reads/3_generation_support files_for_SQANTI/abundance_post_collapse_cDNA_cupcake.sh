#!/bin/bash

# get aboundance of the long reads
singularity exec -B /Beatriz_Toledo/gc /projects/globalscratch/cdna_cupcake.sif get_abundance_post_collapse.py highquality_gmap1_sort_cup.collapsed unpolished.cluster_report.csv
