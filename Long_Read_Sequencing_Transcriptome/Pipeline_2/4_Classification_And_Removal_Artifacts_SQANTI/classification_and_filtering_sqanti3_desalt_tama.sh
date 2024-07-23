#!/bin/bash

#script for classifying transcripts and filtering artifacts using SQANTI3
#SQANTI3 is run as a singularity image file (.sif) within a Singularity container
#use with gtf file generated in step 2 and full-length (full-length_count.txt), isoform (.isoform.results) and splice junction (.twopassSJ.out.tab) quantification obtained in step 3
#in addition, SQANTI3 uses as inputs reference fasta and gtf, cage_peaks and PolyA_motifs files

#load modules
module load apps/singularity/3.5.2

#classify transcripts
singularity exec -B /../SQANTI3_DESALT_TAMA \
    /projects/globalscratch/sqanti3_3.0.sif \
    sqanti3_qc.py /../TAMA/highquality_desalt_sorted_tama.collapsed_nocap.filtered.modified.gtf \
    /../SQANTI3_DESALT_TAMA/Mus_musculus.GRCm38.101.gtf \
    /../SQANTI3_DESALT_TAMA/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    --cage_peak /../SQANTI3_DESALT_TAMA/refTSS_v3.3_mouse_coordinate.mm10.bed \
    --polyA_motif_list /../SQANTI3_DESALT_TAMA/PolyA_motif_List_mm10.txt \ 
    --gtf \
    -e filterdegradationtama_nocap_read_PP_DP_N_PP_rsemSTAR_before_SQANTI_L393.isoforms.results,\
DT.rsemSTAR_L355.isoforms.results,\
DT.rsemSTAR_L222.isoforms.results,\
DT.rsemSTAR_L395.isoforms.results,\
DT.rsemSTAR_L357.isoforms.results,\
DT.rsemSTAR_L224.isoforms.results,\
DT.rsemSTAR_L394.isoforms.results,\
DT.rsemSTAR_L356.isoforms.results,\
DT.rsemSTAR_L223.isoforms.results \
    -fl filterdegradationtama.mapped_fl_count.txt -o PP_DP_N_filtereddegradationama_NOCAP_sqanti2 \
    --isoAnnotLite --gff3 Mus_musculus_GRCm38_Ensembl_86.gff3 \
    -c filterdegradationtamaDP_L356_twopassSJ.out.tab,filterdegradationtamaDP_L223_twopassSJ.out.tab,\
filterdegradationtamaDP_L394_twopassSJ.out.tab,filterdegradationtamaN_L357_twopassSJ.out.tab,\
filterdegradationtamaN_L224_twopassSJ.out.tab,filterdegradationtamaN_L395_twopassSJ.out.tab,\
filterdegradationtamaPP_L355_twopassSJ.out.tab,filterdegradationtamaPP_L222_twopassSJ.out.tab,\
filterdegradationtamaPP_L393_twopassSJ.out.tab

#filter out artifacts and degraded transcripts
singularity exec -B /group/crtd_calegari/People/Beatriz_Toledo/DESALTTIRDTRY/Sqantipreli \
    /projects/globalscratch/sqanti3_3.0.sif \
    python /SQANTI3/sqanti3_RulesFilter.py PP_DP_N_filtereddegradationama_NOCAP_sqanti2_classification.txt \
    PP_DP_N_filtereddegradationama_NOCAP_sqanti2_corrected.fasta \
    PP_DP_N_filtereddegradationama_NOCAP_sqanti2_corrected.gtf
