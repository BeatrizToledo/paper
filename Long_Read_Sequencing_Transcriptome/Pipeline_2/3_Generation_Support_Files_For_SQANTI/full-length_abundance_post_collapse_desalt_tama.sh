#!/bin/bash

#script for generating full-length abundance files for SQANTI3
#use with the collapsed and filtered file highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf from step 2 

#load modules
module load apps/R/4.0.0
module load apps/python/2.7.0 

#get the support levels of the transcripts 
#_trans_read.bed file shows transcript model for each read based on the mapping prior to collapsing. was obtained in the collapsing from step 2
#unpolished.cluster_report.csv of filelist_trans.txt was obtained during Isoseq3 processing

python tama_read_support_levels.py -f /../TAMA/filelist_trans.txt -o /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap -m /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap_trans_read.bed -mt tama

#in R 

Rscript full-length_abundance_post_collapse_desalt_tama.r

#load modules
git clone https://github.com/Magdoll/cDNA_Cupcake.git
export PATH=$PATH:<path_to_Cupcake>/sequence/
export PATH=$PATH:<path_to_Cupcake>/rarefaction/
export PYTHONPATH=$PYTHONPATH:<path_to_Cupcake>/sequence
export PATH=$PATH:<path_to_Cupcake>/sequence
module load apps/python3/3.7.0

#use .fa file from step 2
#scripts from cdna_Cupcake repository 
python /../cDNA_Cupcake/sequence/fa2fq.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.fa
python /../cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq  /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.fastq --read_stat /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered_read_support_modified.txt --classify_csv /../long_reads_fasta/flnc.report.csv -o /../SQANTI_INPUT_DESALT_TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered_full-length_count.txt
