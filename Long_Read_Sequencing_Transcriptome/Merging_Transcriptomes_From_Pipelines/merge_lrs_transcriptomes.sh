#!/bin/bash

#script to merge the transcriptome obtained from pipeline_1 (LRS_1_classification.filtered_lite.gtf) and pipeline_2 (LRS_2_classification.filtered_lite.gtf)
#transcriptomes to be merged are listed in mergelisttama4.txt, result is LRS.gtf

#load modules
module load apps/python/2.7.0 #to check

#merge transcriptome using tama_merge
python tama_merge.py -f mergelistlrs.txt -p LRS.gtf -a 25 -m 10 -z 10 -d merge_dup
