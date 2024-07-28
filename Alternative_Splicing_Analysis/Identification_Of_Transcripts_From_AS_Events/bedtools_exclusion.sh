#!/bin/bash

#script for check which AS events coordinates is located between the "genomic_start_coord" and "genomic_end_coord" from junction file with corresponding isoform

bedtools intersect -a LRS_SRS_alljncmorethan3_coord.txt -b NODES_coord.txt -wo -s -F 0.99 -filenames > LRS_SRS_excl_event_tr.txt
bedtools intersect -a ensembl_alljncmorethan3_coord.txt -b NODES_coord.txt -wo -s -F 0.99 -filenames > ensembl_excl_event_tr.txt
bedtools intersect -a refseq_alljncmorethan3_coord.txt -b NODES_coord.txt -wo -s -F 0.99 -filenames > refseq_excl_event_tr.txt

