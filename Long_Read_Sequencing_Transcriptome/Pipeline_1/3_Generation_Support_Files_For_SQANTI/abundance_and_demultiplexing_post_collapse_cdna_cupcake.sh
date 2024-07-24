#!/bin/bash

#script for generating of demultiplexed files, including abundance files for SQANTI3
#use with the collapsed files  from step 2 
#cDNA Cupcake is run as a singularity image file (.sif) within a Singularity container

#load modules
module load apps/singularity/3.7.1 

#get the fasta file from .gtf file with gffread
gffread -w /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.fasta -g /../CUPCAKE/Mus_musculus.GRCm38.dna.primary_assembly.fa /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.gff

#get abundance of the transcripts
#command generates multiple files, including .read_stat.txt and .abundance.txt (total abundance of transcript)
#unpolished.cluster_report.csv  was obtained during Isoseq3 processing
singularity exec -B /../CUPCAKE /projects/globalscratch/cdna_cupcake.sif get_abundance_post_collapse.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed /../long_reads_fasta/unpolished.cluster_report.csv

#use .fa file
python fa2fq.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.fasta

#demultiplexing. mapped_fl_count.txt file contains full-length read counts for each mapped, unique, isoform for each sample
#flnc.report.csv was obtained during Isoseq3 processing
python demux_isoseq_with_genome.py --mapped_fafq  /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.fastq --read_stat /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.read_stat.txt --classify_csv /../long_reads_fasta/flnc.report.csv -o /../SQANTI_INPUT_GMAP_CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.mapped_fl_count.txt

#generate demultiplexed GFF and FASTA/FASTQ files 
python demux_by_barcode_groups.py /../CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.gff /../SQANTI_INPUT_GMAP_CUPCAKE/highquality_gmap_sorted_cupcake.collapsed.mapped_fl_count.txt /../SQANTI_INPUT_GMAP_CUPCAKE/output_demux_highquality_gmap_sorted_cupcake.collapsed "('bc1001_5p--bc1001_3p','NSC'),('bc1002_5p--bc1002_3p','NP'),('bc1003_5p--bc1003_3p','N’)”



