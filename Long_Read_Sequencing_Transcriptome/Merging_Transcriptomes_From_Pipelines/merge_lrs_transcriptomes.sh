#!/bin/bash

#script to merge the transcriptome obtained from pipeline_1 (LRS_1_classification.filtered_lite.gtf) and pipeline_2 (LRS_2_classification.filtered_lite.gtf)
#mergelist_lrs_1_2.txt contains name of transcriptome files, if samples were capped, their merging priority and the source name

#load modules
module load apps/python/2.7.0
module load apps/R/4.0.0

#make directory for output
mkdir MERGED_TRANSCRIPTOME

#convert .gtf files into .bed files
python /../tama/tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py  /../SQANTI_GMAP_CUPCAKE/LRS_1_classification.filtered_lite.gtf  /../SQANTI_GMAP_CUPCAKE/LRS_1.bed
python /../tama/tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py  /../SQANTI_DESALT_TAMA/LRS_2_classification.filtered_lite.gtf  /../SQANTI_DESALT_TAMA/LRS_2.bed

#merge transcriptome using tama_merge
python /../tama/tama_merge.py -f /../MERGED_TRANSCRIPTOME/mergelist_lrs_1_2.txt -p /../MERGED_TRANSCRIPTOME/LRS_merged -a 25 -m 10 -z 10 -d merge_dup

#convert .bed files into .gtf files
python /../tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py /../MERGED_TRANSCRIPTOME/LRS_merged.bed /../MERGED_TRANSCRIPTOME/LRS_merged.gtf

#run R script to modify .gtf to be compatible with other other tools, such as SQANTI - output LRS.gtf
Rscript /../MERGED_TRANSCRIPTOME/modify_lrs_merged_gtf.r
