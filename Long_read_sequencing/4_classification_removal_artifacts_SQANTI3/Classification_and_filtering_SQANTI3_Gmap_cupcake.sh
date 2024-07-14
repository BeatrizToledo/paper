#!/bin/bash

# Classification_and_filtering_SQANTI3

module load apps/singularity/3.5.2

singularity exec -B  /Beatriz_Toledo/gc/Sqanti /projects/globalscratch/sqanti3_3.0.sif sqanti3_qc.py highquality_gmap1_sort_cup.collapsed.filtered.gff Mus_musculus.GRCm38.101.gtf Mus_musculus.GRCm38.dna.primary_assembly.fa --cage_peak 1refTSS_v3.3_mouse_coordinate.mm10.bed --polyA_motif_list PolyA_motif.txt  --gtf -e GC.rsemSTAR_L393.isoforms.results,GC.rsemSTAR_L355.isoforms.results,GC.rsemSTAR_L222.isoforms.results,GC.rsemSTAR_L395.isoforms.results,GC.rsemSTAR_L357.isoforms.results,GC.rsemSTAR_L224.isoforms.results,GC.rsemSTAR_L394.isoforms.results,GC.rsemSTAR_L356.isoforms.results,GC.rsemSTAR_L223.isoforms.results -fl highquality_gmap1_sort_cup.collapsed.filtered.abundance.txt -o GC_sqanti2  --isoAnnotLite --gff3 Mus_musculus_GRCm38_Ensembl_86.gff3 -c L356_twopassGC.out.tab,L223_twopassGC.out.tab,L394_twopassGC.out.tab,L357_twopassGC.out.tab,L224_twopassGC.out.tab,L395_twopassGC.out.tab,L355_twopassGC.out.tab,L222_twopassGC.out.tab,L393_twopassGC.out.tab

singularity exec -B /Beatriz_Toledo/gc/Sqanti /projects/globalscratch/sqanti3_3.0.sif python /SQANTI3/sqanti3_RulesFilter.py GC_sqanti2_classification.txt GC_sqanti2_corrected.fasta GC_sqanti2_corrected.gtf
