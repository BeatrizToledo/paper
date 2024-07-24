#!/bin/bash

#script for classifying transcripts and filtering artifacts using SQANTI3
#SQANTI3 is run as a singularity image file (.sif) within a Singularity container
#use with gtf file generated in step 2 and full-length (full-length_count.txt), isoform (.isoform.results) and splice junction (.twopassSJ.out.tab) quantification obtained in step 3
#in addition, SQANTI3 uses as inputs reference fasta and gtf, cage_peaks (refTSS database) and PolyA_motifs (supplied in the SQANTI3 data folder) files

#load modules
module load apps/singularity/3.5.2

# #classify transcripts
singularity exec -B /../SQANTI/GMAP_CUPCAKE/ \
    /projects/globalscratch/sqanti3_3.0.sif \
    sqanti3_qc.py /../CUPCAKE/highquality_gmap_sort_cup.collapsed.filtered.gff \
    /../SQANTI_INPUT_GMAP_CUPCAKE/Mus_musculus.GRCm38.101.gtf \
    /../SQANTI_INPUT_GMAP_CUPCAKE/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    --cage_peak /../SQANTI_INPUT_GMAP_CUPCAKE/1refTSS_v3.3_mouse_coordinate.mm10.bed \
    --polyA_motif_list /../SQANTI_INPUT_GMAP_CUPCAKE/PolyA_motif.txt \
    --gtf \
    -e /../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L393.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L355.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L222.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L395.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L357.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L224.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L394.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L356.isoforms.results,\
/../SQANTI_INPUT_GMAP_CUPCAKE/GC.rsemSTAR_L223.isoforms.results \
    -fl /../SQANTI_INPUT_GMAP_CUPCAKE/highquality_gmap_sorted_cupcake.collapsed_full-length_count.txt \
    -o LRS1 \ 
    --isoAnnotLite \
    --gff3 /../SQANTI_INPUT_GMAP_CUPCAKE/Mus_musculus_GRCm38_Ensembl_86.gff3 \
    -c /../SQANTI_INPUT_GMAP_CUPCAKE/L356_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L223_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L394_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L357_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L224_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L395_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L355_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L222_twopassGC.out.tab,\
/../SQANTI_INPUT_GMAP_CUPCAKE/L393_twopassGC.out.tab

#filter out artifacts and degraded transcripts #check folders and final names
singularity exec -B /Beatriz_Toledo/gc/Sqanti \
    /projects/globalscratch/sqanti3_3.0.sif \
    python /SQANTI3/sqanti3_RulesFilter.py GC_sqanti2_classification.txt \
    GC_sqanti2_corrected.fasta \
    GC_sqanti2_corrected.gtf #rename LRS1.gtf?
