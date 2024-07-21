#!/bin/bash

#code to identify and quantify alternative splicing (AS) events using whippet
#uses as input .gtf file derived from the merging of gtf files obtained from long- and short-reads sequencing
#to allow further discovery of cryptic events, an alignment file derived from the merging of single short-read replicate alignments to the merged .gtf file is provided

# load the modules
module load apps/julia/1.6.2

#build an index
julia bin/whippet-index.jl --fasta /../Mus_musculus.GRCm38.dna.primary_assembly.fa --bam /../All.merged.sorted.rmdup.bam --bam-min-reads 10 --gtf /../mergedpcstr2filtereddegradation.nofus_NOCAP_sqanti_corrected_nofus.gtf -x /../Whippet_results

#quantify event percentage splice in (PSI) for each replicate using short-read paired-end sequencing data (_R1.fastq.gz and _R2.fastq.gz). Outputs a .psi.gz file for each replicate)
#sample1 (NSC) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L222_R1.fastq.gz /../fastq/L222_R2.fastq.gz --biascorrect -o /../Whippet_results/L222.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L222.whippet.sam
julia bin/whippet-quant.jl /../fastq/L355_R1.fastq.gz /../fastq/L355_R2.fastq.gz --biascorrect -o /../Whippet_results/L355.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L355.whippet.sam
julia bin/whippet-quant.jl /../fastq/L393_R1.fastq.gz /../fastq/L393_R2.fastq.gz --biascorrect -o /../Whippet_results/L393.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L393.whippet.sam
#sample2 (NP) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L223_R1.fastq.gz /../fastq/L223_R2.fastq.gz --biascorrect -o/../Whippet_results/L223.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L223.whippet.sam
julia bin/whippet-quant.jl /../fastq/L356_R1.fastq.gz /../fastq/L356_R2.fastq.gz --biascorrect -o/../Whippet_results/L356.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L356.whippet.sam
julia bin/whippet-quant.jl /../fastq/L394_R1.fastq.gz /../fastq/L394_R2.fastq.gz --biascorrect -o/../Whippet_results/L394.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L394.whippet.sam
#sample3 (N) - 3 replicates
julia bin/whippet-quant.jl /../fastq/L224_R1.fastq.gz /../fastq/L224_R2.fastq.gz --biascorrect -o/../Whippet_results/L224.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L224.whippet.sam
julia bin/whippet-quant.jl /../fastq/L357_R1.fastq.gz /../fastq/L357_R2.fastq.gz --biascorrect -o/../Whippet_results/L357.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L357.whippet.sam
julia bin/whippet-quant.jl /../fastq/L395_R1.fastq.gz /../fastq/L395_R2.fastq.gz --biascorrect -o/../Whippet_results/L395.whippet -x /../Whippet_results/Whippet_results.jls --sam >/../Whippet_results/L395.whippet.sam

#calculate events differential inclusion (DPSI) between samples
#between NP and NSC
julia bin/whippet-delta.jl -a L223.whippet.psi.gz,L356.whippet.psi.gz, L394.whippet.psi.gz -b L222.whippet.psi.gz,L355.whippet.psi.gz,L393.whippet.psi.gz -s 3 -r 20 -o /../Whippet_results/NP.NSC.DPSI.diff
#between N and NP
julia bin/whippet-delta.jl -a L224.whippet.psi.gz,L357.whippet.psi.gz, L395.whippet.psi.gz -b L223.whippet.psi.gz,L356.whippet.psi.gz,L394.whippet.psi.gz -s 3 -r 20 -o /../Whippet_results/N.NP.DPSI.diff

