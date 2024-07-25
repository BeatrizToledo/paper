#!/bin/bash

#script to merge the transcriptome obtained from pipeline_1 (LRS_1_classification.filtered_lite.gtf) and pipeline_2 (LRS_2_classification.filtered_lite.gtf)
#mergelistlrs.txt contains name of transcriptome files, if samples were capped, their merging priority and the source name

#load modules
module load apps/python/2.7.0

#make directory for output
mkdir MERGED_TRANSCRIPTOME

#convert .gtf files into .bed files
python tama_format_gff_to_bed12_cupcake.py  PP_DP_N_filtereddegradationama_NOCAP_sqanti2_classification.filtered_lite.gtf  PP_DP_N_filtereddegradationama_NOCAP_sqanti2_classification.filtered_lite_nofus.gtf.bed

#merge transcriptome using tama_merge
python /../tama/tama_merge.py -f /../MERGED_TRANSCRIPTOME/mergelistlrs.txt -p /../MERGED_TRANSCRIPTOME/LRS -a 25 -m 10 -z 10 -d merge_dup

python tama_convert_bed_gtf_ensembl_no_cds.py mergedpc2filterdegradationtama.nofus.gtf.bed mergedpc2filtereddegradationama.nofus.gtf
