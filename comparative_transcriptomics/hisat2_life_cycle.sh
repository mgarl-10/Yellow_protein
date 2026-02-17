#!/bin/bash
#$ -S /bin/bash
#$ -N hisat2
#$ -cwd

reads=RNAseq/Trimmed_reads
index=C_alternans

for R1 in ${reads}/*_trimmed_1.fastq.gz; do
    sample=$(basename "$R1" _trimmed_1.fastq.gz)
    R2="${reads}/${sample}_trimmed_2.fastq.gz"

    echo "Running ${sample}"

    hisat2 \
      -x "$index" \
      -1 "$R1" \
      -2 "$R2" \
      -S "${sample}.sam" \
      --rna-strandness RF \
      -p 16
done
