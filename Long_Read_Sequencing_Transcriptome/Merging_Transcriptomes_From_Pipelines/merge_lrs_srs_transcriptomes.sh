#!/bin/bash

#script to merge the transcriptome obtained from LRS (/../MERGED_TRANSCRIPTOME/LRS_merged.bed) and SRS (__)
#mergelistlrssrs.txt contains name of transcriptome files, if samples were capped, their merging priority and the source name

#load modules
module load apps/python/2.7.0

#convert .gtf files into .bed files
python /../tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py  /../stringtie_.gtf /../stringtie__.bed

#merge transcriptome using tama_merge
# give priority to stringtie for the exon boundaries
python /../tama/tama_merge.py -f /../mergelist_lrs_srs.txt -p /../MERGED_TRANSCRIPTOME/LRS_SRS_merged -a 25 -m 10 -z 10 -d merge_dup

#convert .bed files into .gtf files
python /../tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py /../MERGED_TRANSCRIPTOME/LRS_SRS_merged.bed /../MERGED_TRANSCRIPTOME/LRS_SRS_merged.gtf

#run R script to modify .gtf to be compatible with other other tools, such as SQANTI
Rscript /../MERGED_TRANSCRIPTOME/modify_lrs_srs_merged_gtf.r
