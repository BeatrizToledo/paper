#!/bin/bash

#In order to quantify with whippet I need to first generate a single BAM file from the differen cell types
#prepare bam for Whippet
#load the modules

qlogin -q all.q -l h_rt=44:00:00 -l mem_free=70G -R yes -now n -pe smp 1


module load apps/samtools/1.9

samtools merge -@ 70 hisat.bam L222.bam L223.bam L224.bam L355.bam L356.bam L357.bam L393.bam L394.bam L395.bam
samtools sort -@ 70 -o hisat.sort hisat.bam
samtools rmdup -S hisat.sort hisat.sort.rmdup.bam
samtools index hisat.sort.rmdup.bam


#   Whippet was developed with Julia

module load apps/julia/1.6.2 
cd Whippet.jl/
julia --project -e 'using Pkg; Pkg.instantiate()'


# Generation of index
julia bin/whippet-index.jl --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa --bam hisat.sort.rmdup.bam --bam-min-reads 10 --gtf PCST.gtf -x PCSTindex 

# Quantification in each individual cell 
julia bin/whippet-quant.jl L222_R1.fastq L222_R2.fastq --biascorrect -x  PCSTindex.jls -o L222
julia bin/whippet-quant.jl L223_R1.fastq L223_R2.fastq --biascorrect -x  PCSTindex.jls -o L223
julia bin/whippet-quant.jl L224_R1.fastq L224_R2.fastq --biascorrect -x  PCSTindex.jls -o L224
julia bin/whippet-quant.jl L394_R1.fastq L394_R2.fastq --biascorrect -x  PCSTindex.jls -o L394
julia bin/whippet-quant.jl L393_R1.fastq L393_R2.fastq --biascorrect -x  PCSTindex.jls -o L393
julia bin/whippet-quant.jl L356_R1.fastq L356_R2.fastq --biascorrect -x  PCSTindex.jls -o L356
julia bin/whippet-quant.jl L355_R1.fastq L355_R2.fastq --biascorrect -x  PCSTindex.jls -o L355
julia bin/whippet-quant.jl L357_R2.fastq L357_R1.fastq --biascorrect -x  PCSTindex.jls -o L357
julia bin/whippet-quant.jl L395_R1.fastq L395_R2.fastq --biascorrect -x  PCSTindex.jls -o L395 


# calculation of differential alternative splicing
julia bin/whippet-delta.jl -a L223.psi.gz,L356.psi.gz,L394.psi.gz -b L222.psi.gz,L355.psi.gz,L393.psi.gz -s 3 -o PP.DP.DPSI.diff
julia bin/whippet-delta.jl -a L224.psi.gz,L357.psi.gz,L395.psi.gz -b L223.psi.gz,L356.psi.gz,L394.psi.gz -s 3 -o DP.N.DPSI.diff




