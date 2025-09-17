#!/bin/bash
#$ -S /bin/bash
#$ -N featurec
#$ -cwd

/ebio/ag-salem/projects/CassidinaeGenomics/code/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 16 -s 2 -t CDS,rRNA,misc_RNA,tRNA -g Name -a Stammera_Chelymorpha_alternans.gff -o Stammera_featurecounts_all_samples_single_Name.txt CTL_UNT_R1_Stammera_single.sam CTL_HUM_R1_Stammera_single.sam CTL_UNT_R2_Stammera_single.sam CTL_HUM_R2_Stammera_single.sam CTL_UNT_R3_Stammera_single.sam CTL_HUM_R3_Stammera_single.sam GFP_UNT_R1_Stammera_single.sam GFP_HUM_R1_Stammera_single.sam GFP_UNT_R2_Stammera_single.sam GFP_HUM_R2_Stammera_single.sam GFP_UNT_R3_Stammera_single.sam GFP_HUM_R3_Stammera_single.sam YELL_UNT_R1_Stammera_single.sam YELL_HUM_R1_Stammera_single.sam YELL_UNT_R2_Stammera_single.sam YELL_HUM_R2_Stammera_single.sam YELL_UNT_R3_Stammera_single.sam YELL_HUM_R3_Stammera_single.sam