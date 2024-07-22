#!/bin/bash

#use with the alignment file from deSALT highquality_desalt_sorted.sam (from step 1)

#load modules
module load apps/python/2.7.0 
module load apps/R/4.0.0
module load apps/gffread/0.12.1
git clone https://github.com/GenomeRIK/tama

#change chromosome nomenclature removing "chr" to make the .sam file compatible with tama
cat /../DESALT/highquality_desalt.sam | sed 's/chrM/MT/g' > /../DESALT/highquality_desalt_sorted.sam
cat /../DESALT/highquality_desalt_sorted.sam | sed 's/chr//g' > /../DESALT/highquality_desalt_sorted.sam

#make directory for tama outputs
mkdir TAMA

#Collapsing step. tama_collapse.py command generates 9 files, including .bed file.
#samples were not capped, still generate file with -x capped and -x no_cap to evaluate the degration level. 
python /../tama/tama_collapse.py -s /../DESALT/highquality_desalt_sorted.sam -f Mus_musculus.GRCm38.dna.primary_assembly.fa -p  /../TAMA/highquality_desalt_sorted_tama.collapsed_capped -x capped 
python /../tama/tama_collapse.py -s /../DESALT/highquality_desalt_sorted.sam -f Mus_musculus.GRCm38.dna.primary_assembly.fa -p  /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap -x no_cap

#calculate the degradation signature with .bed files generated from tama_collapse.py command
python /../tama/tama_go/file_stats/tama_degradation_signature.py -c /../TAMA/highquality_desalt_sorted_tama.collapsed_capped.bed -nc /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.bed  -o /../TAMA/highquality_desalt_tama_collapsed_DegSig

#filter out degraded transcripts with a max difference of 100 nt in the 3' end (option -e 100); keep the longer transcript model 
python /../tama/tama_go/filter_transcript_models/tama_remove_fragment_models.py -f /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.bed -o /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.bed -e 100

#get gtf file from the bed file
python /../tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.bed /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.gtf

#!!!IMPORTANT STEP!!!
# The gtf file generated by tama has a different structure than the one generated by cdna_cupcake
# tama gtf 
#SQANTI accepts gtf with formatting as the one of cupcake. The gtf of tama must be modified as follows to make it compatible with SQANTI
# 1       PBRI  gene      3199736 3202443 .       -       .       gene_id "G1.1";
# 1       PBRI  transcript        3199736 3202443 .       -       .       gene_id "G1.1"; transcript_id "transcript/132038"; uniq_trans_id "transcript/132038";
# 1       PBRI  exon      3199736 3202443 .       -       .       gene_id "G1.1"; transcript_id "transcript/132038"; exon_number "1"; uniq_trans_id "transcript/132038";
# cupcake gtf 
#1       PacBio  transcript      3136606 3136668 .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
#1       PacBio  exon    3136606 3136668 .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";

#run R script to modify table
Rscript Modify_Tama_table.r

#get the fasta file from the modified gtf file with gffread
gffread -w /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.fa -g /../TAMA/Mus_musculus.GRCm38.dna.primary_assembly.fa /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf
