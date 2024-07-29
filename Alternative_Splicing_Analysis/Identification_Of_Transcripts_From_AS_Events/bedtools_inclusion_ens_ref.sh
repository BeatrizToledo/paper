#!/bin/bash

#script for checking which of the remaining AS events coordinates overlaps an exon

bedtools intersect -a /../ISOFORM/ensembl_gtf_alljncmorethan3_coord.txt -b /../ISOFORM/as_not_overlap_coord.txt -wo -s -F 0.99 -filenames > as_overlap_ens.txt
bedtools intersect -a /../ISOFORM/refseq_gtf_alljncmorethan3_coord.txt -b /../ISOFORM/as_not_overlap_coord.txt -wo -s -F 0.99 -filenames > as_overlap_ref.txt
