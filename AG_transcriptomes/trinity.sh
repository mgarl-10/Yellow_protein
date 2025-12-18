#!/bin/bash
#$ -S /bin/bash
#$ -N trinity_Asp
#$ -cwd

source ~/.bashrc 

forward=(*trimmed_1.fastq.gz)
reverse=(*trimmed_2.fastq.gz)

for i in ${!forward[@]};do
  f=${forward[i]}
  r=${reverse[i]}
  tag=${f%*trimmed_1.fastq.gz}
Trinity --seqType fq --CPU 64 --max_memory=500G --normalize_by_read_set --output "$tag"_trinity_out --left "$f" --right "$r" --SS_lib_type RF
done
