#!/bin/bash

#script for classifying transcripts and filtering fusion transcripts using SQANTI3
#SQANTI3 is run as a singularity image file (.sif) within a Singularity container
#use final merged transcriptome (LRS_SRS.gtf)

#load modules
module load apps/singularity/3.7.1
module load apps/R/4.0.0

#classification of transcripts from LRS_SRS
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py \
/../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.101.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--cage_peak /../REFERENCE_INPUT/1refTSS_v3.3_mouse_coordinate.mm10.bed --polyA_motif_list  /../REFERENCE_INPUT/PolyA_motif_List_mm10.txt  \
--gtf -o /../MERGED_TRANSCRIPTOME/LRS_SRS  --isoAnnotLite --gff3 /../REFERENCE_INPUT/Mus_musculus_GRCm38_Ensembl_86.gff3 

#classification with refseq and gencode
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/mm10.ncbiRefSeq2.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o LRS_SRS_refseq --skipORF --skip_report
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/gencode.vM10.primary_assembly.annotation.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o LRS_SRS_gencode --skipORF —skip_report

#remove fusion transcripts from LRS_SRS
Rscript /../MERGED_TRANSCRIPTOME/classification_transcripts_and_filtering_of_fusion_transcripts.r

