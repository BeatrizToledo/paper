#!/bin/bash

#code to identify and quantify alternative splicing (AS) events using whippet
#uses as input .gtf file derived from the merging of gtf files obtained from long- and short-reads sequencing
#to allow further discovery of cryptic events, an alignment file derived from the merging of single short-read replicate alignments to the merged .gtf file is provided

# load the modules
module load apps/julia/1.6.2

#build an index
julia bin/whippet-index.jl --fasta /../Mus_musculus.GRCm38.dna.primary_assembly.fa --bam /../All.merged.sorted.rmdup.bam --bam-min-reads 10 --gtf /../mergedpcstr2filtereddegradation.nofus_NOCAP_sqanti_corrected_nofus.gtf -x /../Whippet_results

#quantify event percentage splice in (PSI) for each replicate using short-read sequencing data (output .psi.gz file for each sample)
#sample1 (NSC) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L222_rfp_minus_Track_3438_3447_R1.fastq.gz /../fastq/L222_rfp_minus_Track-3438-3447_R2.fastq.gz --biascorrect -o /../Whippet_results/L222.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L222.whippet.sam
julia bin/whippet-quant.jl /../fastq/L355_rfp_minus_Track_3441_3450_R1.fastq.gz /../fastq/L355_rfp_minus_Track-3441-3450_R2.fastq.gz --biascorrect -o /../Whippet_results/L355.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L355.whippet.sam
julia bin/whippet-quant.jl /../fastq/L393_rfp_minus_Track_3444_3453_R1.fastq.gz /../fastq/L393_rfp_minus_Track-3444-3453_R2.fastq.gz --biascorrect -o /../Whippet_results/L393.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L393.whippet.sam
#sample2 (NP) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L223_rfp_plus_Track_3439_3448_R1.fastq.gz /../fastq/L223_rfp_plus_Track-3439-3448_R2.fastq.gz --biascorrect -o/../Whippet_results/L223.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L223.whippet.sam
julia bin/whippet-quant.jl /../fastq/L356_rfp_plus_Track_3442_3451_R1.fastq.gz /../fastq/L356_rfp_plus_Track-3442-3451_R2.fastq.gz --biascorrect -o/../Whippet_results/L356.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L356.whippet.sam
julia bin/whippet-quant.jl /../fastq/L394_rfp_plus_Track_3445_3454_R1.fastq.gz /../fastq/L394_rfp_plus_Track-3445-3454_R2.fastq.gz --biascorrect -o/../Whippet_results/L394.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L394.whippet.sam
#sample3 (N) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L224_gfp_plus_Track_3440_3449_R1.fastq.gz /../fastq/L224_gfp_plus_Track-3440-3449_R2.fastq.gz --biascorrect -o/../Whippet_results/L224.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L224.whippet.sam
julia bin/whippet-quant.jl /../fastq/L357_gfp_plus_Track_3443_3452_R1.fastq.gz /../fastq/L357_gfp_plus_Track-3443-3452_R2.fastq.gz --biascorrect -o/../Whippet_results/L357.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L357.whippet.sam
julia bin/whippet-quant.jl /../fastq/L395_gfp_plus_Track_3446_3455_R1.fastq.gz /../fastq/L395_gfp_plus_Track-3446-3455_R2.fastq.gz --biascorrect -o/../Whippet_results/L395.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L395.whippet.sam

#calculate events differential inclusion (DPSI) between samples
#between NP and NSC
julia bin/whippet-delta.jl -a L223.whippet.psi.gz,L356.whippet.psi.gz, L394.whippet.psi.gz -b L222.whippet.psi.gz,L355.whippet.psi.gz,L393.whippet.psi.gz -s 3 -r 20 -o /../Whippet_results/NP.NSC.DPSI.diff
#between N and NP
julia bin/whippet-delta.jl -a L224.whippet.psi.gz,L357.whippet.psi.gz, L395.whippet.psi.gz -b L223.whippet.psi.gz,L356.whippet.psi.gz,L394.whippet.psi.gz -s 3 -r 20 -o /../Whippet_results/NP.NSC.DPSI.diff

