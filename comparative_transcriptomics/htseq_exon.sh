#!/bin/bash
#$ -S /bin/bash
#$ -N htseq_all
#$ -cwd

source ~/.bashrc

gff=Chelymorpha_alternans.gff3

for bam in ../*_sorted.bam; do
    sample=$(basename "$bam" _sorted.bam)
    out="count_table_${sample}_exon_Parent_new_ann.txt"

    echo "Processing $sample"

    htseq-count \
      -f bam \
      -t exon \
      --stranded=reverse \
      -r pos \
      -i Parent \
      "$bam" \
      "$gff" \
      > "$out"
done
