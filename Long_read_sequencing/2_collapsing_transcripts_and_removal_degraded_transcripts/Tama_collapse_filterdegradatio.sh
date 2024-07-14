#!/bin/bash

#In order to collapse redundant transcripts with tama is first needed to modify the .sam file. Need to remove chr and replace chrM with MT

cat hq.rn.sorted.sam | sed 's/chrM/MT/g' > hq.rn1.sorted.sam
cat hq.rn1.sorted.sam | sed 's/chr//g' > hq.rn2.sorted.sam


# My samples were not capped, still I generated file with -x capped and -x no_cap to evaluate the degration level 

python tama_collapse.py -s hq.rn2.sorted.sam -f Mus_musculus.GRCm38.dna.primary_assembly.fa -p  hq.rn.collapse_capped -x capped 
python tama_collapse.py -s hq.rn2.sorted.sam -f Mus_musculus.GRCm38.dna.primary_assembly.fa -p  hq.rn.collapse_nocap -x no_cap


# To calculate the degradation signature

python tama_degradation_signature.py -c hq.rn.collapse_capped.bed -nc hq.rn.collapse_nocap.bed  -o highquality_desalt_collapse_DegSig


# for filtering degraded transcripts I chose the option -e 100 

tama_remove_fragment_models.py -f filterdegradationtama.bed -o filterdegradationtama1 -e 100

# In order to make the tama gtf look like the cupcake gtf (which is needed for using it as input for SQANTI) I did the following
# tama gtf looks like
# 1       PBRI  gene      3199736 3202443 .       -       .       gene_id "G1.1";
# 1       PBRI  transcript        3199736 3202443 .       -       .       gene_id "G1.1"; transcript_id "transcript/132038"; uniq_trans_id "transcript/132038";
# 1       PBRI  exon      3199736 3202443 .       -       .       gene_id "G1.1"; transcript_id "transcript/132038"; exon_number "1"; uniq_trans_id "transcript/132038";
# Cupcake gtf look like
#1       PacBio  transcript      3136606 3136668 .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
#1       PacBio  exon    3136606 3136668 .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
# First i need to get gtf file from the bed file

python tama_convert_bed_gtf_ensembl_no_cds.py filterdegradationtama.bed filterdegradationtama.gtf

#Run R script to modify table
Rscript Modify_Tama_table.r

#get the fasta file from the modified gtf file with gffread

gffread -w filterdegradationtama.modified.fa -g Mus_musculus.GRCm38.dna.primary_assembly.fa filterdegradationtama.modified.gtf
