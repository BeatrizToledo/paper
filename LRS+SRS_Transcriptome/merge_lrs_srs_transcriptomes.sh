#!/bin/bash

#script to merge the transcriptome obtained from LRS (/../MERGED_TRANSCRIPTOME/LRS_merged.bed) and SRS (SRS.gtf)
#mergelist_lrs_srs.txt contains name of transcriptome files, if samples were capped, their merging priority and the source name

#load modules
module load apps/python/2.7.0
module load apps/R/4.0.0

#convert .gtf files into .bed files
python /../tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py /../ASSEMBLY/SRS.gtf /../ASSEMBLY/SRS.bed

#merge transcriptome using tama_merge
# give priority to stringtie for the exon boundaries
python /../tama/tama_merge.py -f /../mergelist_lrs_srs.txt -p /../MERGED_TRANSCRIPTOME/LRS_SRS -a 25 -m 10 -z 10 -d merge_dup

#convert .bed files into .gtf files
python /../tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py /../MERGED_TRANSCRIPTOME/LRS_SRS.bed /../MERGED_TRANSCRIPTOME/LRS_SRS.gtf

#run R script to modify .gtf to be compatible with other other tools, such as SQANTI - OUTPUT LRS_SRS.gtf
Rscript /../MERGED_TRANSCRIPTOME/modify_lrs_srs_merged_gtf.r
