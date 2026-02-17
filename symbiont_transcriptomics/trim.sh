#!/bin/bash
#$ -S /bin/bash
#$ -N trim
#$ -cwd


forward=(*_R1_001.fastq.gz)
reverse=(*_R2_001.fastq.gz)

for i in ${!forward[@]};do
  f=${forward[i]}
  r=${reverse[i]}
  tag=${f%*_R1_001.fastq.gz}
java -jar /ebio/ag-salem/projects/CassidinaeGenomics/code/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 16 "$f" "$tag"_trimmed_single_1.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
done
