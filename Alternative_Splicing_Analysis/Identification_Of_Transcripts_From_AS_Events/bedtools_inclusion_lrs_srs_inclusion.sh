#!/bin/bash

#script for checking which of the remaining AS events coordinates overlaps an exon

bedtools intersect -a /../ISOFORM/LRS_SRS_alljncmorethan3_coord.txt -b /../ISOFORM/as_not_start_and_or_end_coord.txt -wo -s -F 0.99 -filenames > /../ISOFORM/as_overlap.txt
