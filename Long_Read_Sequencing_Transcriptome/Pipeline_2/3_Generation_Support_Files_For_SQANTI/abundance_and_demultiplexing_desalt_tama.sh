#!/bin/bash

#script for generating of demultiplexed files, including full-length abundance files for SQANTI3
#use with the collapsed file highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf from step 2 

#load modules
module load apps/R/4.0.0
module load apps/python/2.7.0 
module load apps/gffread/0.12.1

#get the fasta file from the gtf file with gffread
gffread -w /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.fa -g /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.gtf

#get the abundance levels of the transcripts 
#_trans_read.bed file shows transcript model for each read based on the mapping prior to collapsing. was obtained in the collapsing from step 2
#unpolished.cluster_report.csv of filelist_trans.txt was obtained during Isoseq3 processing
python /../tama/tama_go/read_support/tama_read_support_levels.py -f /../TAMA/filelist_trans.txt -o /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap -m /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap_trans_read.bed -mt tama

#run R script to modify table
Rscript /../TAMA/full-length_abundance_post_collapse_desalt_tama.r

#load modules
#make sure you have a git clone of cDNA_Cupcake
module unload apps/python3/2.7.0
module load apps/python3/3.7.0

#use .fa file 
#scripts from cDNA_Cupcake repository 
python /../cDNA_Cupcake/sequence/fa2fq.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.fa

#demultiplexing. mapped_fl_count.txt file contains full-length read counts for each mapped, unique, isoform for each sample
#flnc.report.csv was obtained during Isoseq3 processing
python /../cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq  /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.fastq --read_stat /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap_read_support_modified.txt --classify_csv /../long_reads_fasta/flnc.report.csv -o /../SQANTI_INPUT_DESALT_TAMA/highquality_desalt_sorted_tama.collapsed_nocap_full-length_count.txt

#generate demultiplexed GFF and FASTA/FASTQ files 
python /../cDNA_Cupcake/post_isoseq_cluster/demux_by_barcode_groups.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.gff /../SQANTI_INPUT_DESALT_TAMA/highquality_desalt_sorted_tama.collapsed_nocap_full-length_count.txt /../SQANTI_INPUT_DESALT_TAMA/output_demux_highquality_desalt_sorted_tama.collapsed_nocap "('bc1001_5p--bc1001_3p','NSC'),('bc1002_5p--bc1002_3p','NP'),('bc1003_5p--bc1003_3p','N’)”
