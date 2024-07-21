#!/bin/bash

#code to merge the transcriptome obtained from pipeline_1 (LRS_1.gtf) and pipeline_2 (LRS_2.gtf)
#transcriptomes to be merged are listed in mergelisttama4.txt, result is LRS.gtf

#load modules
module load apps/python/2.9 #ro check

#merge transcriptome using tama_merge
python tama_merge.py -f mergelisttama4.txt -p LRS.gtf -a 25 -m 10 -z 10 -d merge_dup
