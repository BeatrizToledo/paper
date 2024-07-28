#!/bin/bash

#script for classifying transcripts and filtering fusion transcripts using SQANTI3
#SQANTI3 is run as a singularity image file (.sif) within a Singularity container
#use final merged transcriptome (LRS_SRS.gtf)

#load modules
module load apps/singularity/3.7.1
module load apps/R/4.0.0

#classification of transcripts from LRS_SRS with ensembl
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py \
/../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.101.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--cage_peak /../REFERENCE_INPUT/1refTSS_v3.3_mouse_coordinate.mm10.bed --polyA_motif_list  /../REFERENCE_INPUT/PolyA_motif_List_mm10.txt  \
--gtf -o /../MERGED_TRANSCRIPTOME/LRS_SRS  --isoAnnotLite --gff3 /../REFERENCE_INPUT/Mus_musculus_GRCm38_Ensembl_86.gff3 -c /../MERGED_TRANSCRIPTOME/L356_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L223_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L394_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L357_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L224_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L395_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L355_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L222_lrs_srs_sj.out.tab,\
/../MERGED_TRANSCRIPTOME/L393_lrs_srs_sj.out.tab

#classification of transcripts from LRS_SRS with refseq and gencode
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/mm10.ncbiRefSeq2.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o /../MERGED_TRANSCRIPTOME/LRS_SRS_refseq --skipORF --skip_report
singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../MERGED_TRANSCRIPTOME/LRS_SRS.gtf /../REFERENCE_INPUT/gencode.vM10.primary_assembly.annotation.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o /../MERGED_TRANSCRIPTOME/LRS_SRS_gencode --skipORF —skip_report

#remove fusion transcripts from LRS_SRS
Rscript /../MERGED_TRANSCRIPTOME/classification_transcripts_and_filtering_of_fusion_transcripts.r

#classification of ensembl, refseq and gencode transcripts

singularity exec -B /../MERGED_TRANSCRIPTOME,/../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../REFERENCE_INPUT/mm10.ncbiRefSeq2.gtf /../REFERENCE_INPUT/mm10.ncbiRefSeq2.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o /../REFERENCE_INPUT/refseq --skipORF --skip_report -c /../REFERENCE_INPUT/L356_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L223_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L394_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L357_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L224_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L395_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L355_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L222_refseq_sj.out.tab,\
/../REFERENCE_INPUT/L393_refseq_sj.out.tab
singularity exec -B /../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../REFERENCE_INPUT/gencode.vM10.primary_assembly.annotation.gtf /../REFERENCE_INPUT/gencode.vM10.primary_assembly.annotation.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o /../REFERENCE_INPUT/gencode --skipORF —skip_report  -c /../REFERENCE_INPUT/L356_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L223_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L394_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L357_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L224_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L395_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L355_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L222_gencode_sj.out.tab,\
/../REFERENCE_INPUT/L393_gencode_sj.out.tab
singularity exec -B /../REFERENCE_INPUT /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py  /../REFERENCE_INPUT/Mus_musculus.GRCm38.101.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.101.gtf /../REFERENCE_INPUT/Mus_musculus.GRCm38.dna.primary_assembly.fa   --gtf -o /../REFERENCE_INPUT/ensembl --skipORF —skip_report  -c /../REFERENCE_INPUT/L356_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L223_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L394_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L357_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L224_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L395_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L355_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L222_ensembl_sj.out.tab,\
/../REFERENCE_INPUT/L393_ensembl_sj.out.tab


