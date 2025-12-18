#!/bin/bash
#$ -S /bin/bash
#$ -N map_sing
#$ -cwd

source ~/.bashrc
conda activate bowtie2

INDEX=Stammera

for read in *_trimmed_single_1.fastq.gz; do
    sample=$(basename "$read" _trimmed_single_1.fastq.gz)

    echo "Mapping $sample"

    bowtie2 \
      -x "$INDEX" \
      -U "$read" \
      --fr \
      --no-unal \
      -p 64 \
      -S "${sample}_Stammera_single.sam"
done
