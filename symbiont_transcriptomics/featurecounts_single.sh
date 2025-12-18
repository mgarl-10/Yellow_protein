#!/bin/bash
#$ -S /bin/bash
#$ -N featurec
#$ -cwd

FEATURECOUNTS=/ebio/ag-salem/projects/CassidinaeGenomics/code/subread-2.0.3-Linux-x86_64/bin/featureCounts

$FEATURECOUNTS \
  -T 16 \
  -s 2 \
  -t CDS,rRNA,misc_RNA,tRNA \
  -g Name \
  -a Stammera_Chelymorpha_alternans.gff \
  -o Stammera_featurecounts_all_samples_single_Name.txt \
  *_Stammera_single.sam
