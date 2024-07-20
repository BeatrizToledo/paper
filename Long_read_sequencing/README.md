Processing of long-read sequencing data with different bioinformatic tools

Includes the steps: 

1_alignement_HQ_transcripts   
2_collapsing_transcripts_and_removal_degraded_transcripts  
3_generation_support files_for_SQANTI  
4_classification_removal_artifacts_SQANTI3  
5_Merging_transcriptomes  
6_Identification_isoform_from_AS  

For step 1) and 2) different tools were used, resulting in two pipelines:  
a) DeSALT (step 1) + Tama (step 2)  
b) GMAP (step 1) + Cupcake (step2)  

The results of the two pipelines were processed separately (steps 3-4) and finally merged (step 5).
