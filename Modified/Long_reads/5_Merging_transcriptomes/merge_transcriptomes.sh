#!/bin/bash

#code to merge the transcriptome obtained from gmap+cdna_cupcake (pipeline 1) or desalt+tama (pipeline 2)
#transcriptomes to be merged are listed in mergelisttama4.txt, result is mergedpc2filterdegradationtama.nofus.gtf

#load modules
module load apps/python

#merge transcriptome using tama_merge
python tama_merge.py -f mergelisttama4.txt -p mergedpc2filterdegradationtama.nofus.gtf -a 25 -m 10 -z 10 -d merge_dup
