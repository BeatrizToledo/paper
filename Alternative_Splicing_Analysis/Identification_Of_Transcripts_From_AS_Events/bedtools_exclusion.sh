#!/bin/bash

#script for check which AS events coordinates is located between the "genomic_start_coord" and "genomic_end_coord" from junction file with corresponding isoform

bedtools intersect -a /../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt -b /../ISOFORM/NODES_coord.txt -wo -s -F 0.99 -filenames > /../ISOFORM/LRS_SRS_excl_event_tr.txt
bedtools intersect -a /../ISOFORM/ensembl_alljncmorethan3_coord.txt -b /../ISOFORM/NODES_coord.txt -wo -s -F 0.99 -filenames > /../ISOFORM/ensembl_excl_event_tr.txt
bedtools intersect -a /../ISOFORM/refseq_alljncmorethan3_coord.txt -b /../ISOFORM/NODES_coord.txt -wo -s -F 0.99 -filenames > /../ISOFORM/refseq_excl_event_tr.txt

