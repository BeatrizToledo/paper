#!/bin/bash

#script for generating full-length abundance files of demultiplexed samples for SQANTI3
#use with the collapsed and filtered files  from step 2 
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 

#get abundance of the transcripts
#command generates multiple files, including .read_stat.txt and .abundance.txt (total abundance of transcript)
singularity exec -B /../CUPCAKE /projects/globalscratch/cdna_cupcake.sif get_abundance_post_collapse.py /../CUPCAKE/highquality_gmap_sorted_cup.collapsed.filtered /../long_reads_fasta/unpolished.cluster_report.csv

#use .fa file from step 2
python fa2fq.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.fasta

#demultiplexing. mapped_fl_count.txt file contains full-length read counts for each mapped, unique, isoform for each sample
python demux_isoseq_with_genome.py --mapped_fafq  /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.fastq --read_stat /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.read_stat.txt --classify_csv /../long_reads_fasta/flnc.report.csv -o /../SQANTI_INPUT_GMAP_CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.mapped_fl_count.txt

#generate demultiplexed GFF and FASTA/FASTQ files 
python demux_by_barcode_groups.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.gff /../SQANTI_INPUT_GMAP_CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.filtered.mapped_fl_count.txt /../SQANTI_INPUT_GMAP_CUPCAKE/output_demux_highquality_gmap_sorted_cupcake.collapsed.filtered "('bc1001_5p--bc1001_3p','PP'),('bc1002_5p--bc1002_3p','DP'),('bc1003_5p--bc1003_3p','N’)”



