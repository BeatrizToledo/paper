#!/bin/bash

#script for classifying transcripts and filtering artifacts using SQANTI3
#SQANTI3 is run as a singularity image file (.sif) within a Singularity container
#use with gtf file generated in step 2 and full-length (full-length_count.txt), isoform (.isoform.results) and splice junction (.twopassSJ.out.tab) quantification obtained in step 3
#in addition, SQANTI3 uses as inputs reference fasta and gtf, cage_peaks (refTSS database) and PolyA_motifs (supplied in the SQANTI3 data folder) files

#load modules
module load apps/singularity/3.7.1

#classify transcripts
singularity exec -B /../SQANTI3_DESALT_TAMA,/../TAMA,/../SQANTI_INPUT_DESALT_TAMA \
    /projects/globalscratch/sqanti3_3.0.sif \
    sqanti3_qc.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf \
    /../SQANTI_INPUT_DESALT_TAMA/Mus_musculus.GRCm38.101.gtf \
    /../SQANTI_INPUT_DESALT_TAMA/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    --cage_peak /../SQANTI_INPUT_DESALT_TAMA/refTSS_v3.3_mouse_coordinate.mm10.bed \
    --polyA_motif_list /../SQANTI_INPUT_DESALT_TAMA/PolyA_motif_List_mm10.txt \ 
    --gtf \
    -e /../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L393.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L355.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L222.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L395.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L357.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L224.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L394.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L356.isoforms.results,\
/../SQANTI_INPUT_DESALT_TAMA/dt.rsemstar_L223.isoforms.results \
    -fl /../SQANTI_INPUT_DESALT_TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered_full-length_count.txt -o /../SQANTI_DESALT_TAMA/LRS2 \
    --isoAnnotLite --gff3 /../SQANTI_INPUT_DESALT_TAMA/Mus_musculus_GRCm38_Ensembl_86.gff3 \
    -c /../SQANTI_INPUT_DESALT_TAMA/L356_dt_sj.out.tab,/../SQANTI_INPUT_DESALT_TAMA/L223_dt_sj.out.tab,\
/../SQANTI_INPUT_DESALT_TAMA/L394_dt_sj.out.tab,/../SQANTI_INPUT_DESALT_TAMA/L357_dt_sj.out.tab,\
/../SQANTI_INPUT_DESALT_TAMA/L224_dt_sj.out.tab,/../SQANTI_INPUT_DESALT_TAMA/L395_dt_sj.out.tab,\
/../SQANTI_INPUT_DESALT_TAMA/L355_dt_sj.out.tab,/../SQANTI_INPUT_DESALT_TAMA/L222_dt_sj.out.tab,\
/../SQANTI_INPUT_DESALT_TAMA/L393_dt_sj.out.tab

#filter out artifacts and degraded transcripts
singularity exec -B /../SQANTI3_DESALT_TAMA \
    /projects/globalscratch/sqanti3_3.0.sif \
    python /../SQANTI_DESALT_TAMA/sqanti3_RulesFilter.py /SQANTI_DESALT_TAMA/LRS2_classification.txt \
    /../SQANTI_DESALT_TAMA/LRS2_corrected.fasta \
    /../SQANTI_DESALT_TAMA/LRS2_corrected.gtf
