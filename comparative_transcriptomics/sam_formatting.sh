#!/bin/bash
#$ -S /bin/bash
#$ -N samtools
#$ -cwd

for sam in *.sam; do
    sample=$(basename "$sam" .sam)

    echo "Processing $sample"

    samtools view -b "$sam" | \
    samtools sort -o "${sample}_sorted.bam"
done
