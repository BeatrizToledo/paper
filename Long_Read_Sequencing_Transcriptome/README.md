This folder cointains the codes for processing of long-read RNA sequencing data to obtain sample-specific transcriptomes. 

Includes the steps: 

1_Alignement  
2_Collapsing_Transcripts_And_Removal_Degraded_Transcripts  
3_Generation_Support_Files_For_SQANTI  
4_Classification_And_Removal_Artifacts_SQANTI  

For step 1) and 2) different tools were used, resulting in two pipelines:  

Pipeline_1: GMAP (step 1) + cDNA_Cupcake (step2) 
Pipeline_2: deSALT (step 1) + Tama (step 2)  

The results of the two pipelines were processed separately (steps 3-4) and finally merged (see Merging_Transcriptomes_From_Pipelines) into a single long_read_sequencing_transcriptome (LRS.gtf).

The data come from bulk long-read sequencing of sorted neural stem cells (NSC), committed neural progenitors (NP) and newborn neurons (N) from lateral neocortices of developing mouse brains at embryonic day (E) 14.5.
