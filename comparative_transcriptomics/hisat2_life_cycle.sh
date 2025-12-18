#!/bin/bash
#$ -S /bin/bash
#$ -N hisat2_LC
#$ -cwd

HISAT2=/ebio/ag-salem/projects/CassidinaeGenomics/code/hisat2-2.2.0/hisat2
READS=/ebio/ag-salem/projects/BeetleGenomes/data/RNAseq/Stammera_RNAseq/Trimmed_reads
INDEX=C_alternans

for i in {1..12}; do
    echo "Running sample S${i}"

    $HISAT2 \
      -x $INDEX \
      -1 ${READS}/S${i}_trimmed_1.fastq.gz \
      -2 ${READS}/S${i}_trimmed_2.fastq.gz \
      -S S${i}.sam \
      --rna-strandness RF \
      -p 16
done
