#!/bin/bash
#$ -S /bin/bash
#$ -N htseq_lifecyle
#$ -cwd

source ~/.bashrc


htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_FG_1_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_FG_1_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_FG_2_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_FG_2_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_FG_3_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_FG_3_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_OV_1_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_OV_1_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_OV_2_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_OV_2_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Adult_OV_3_sorted.bam Chelymorpha_alternans.gff3 > count_table_Adult_OV_3_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Larvae_1_sorted.bam Chelymorpha_alternans.gff3 > count_table_Larvae_1_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Larvae_2_sorted.bam Chelymorpha_alternans.gff3 > count_table_Larvae_2_exon_Parent_new_ann.txt

htseq-count -f bam -t exon --stranded=reverse -r pos -i Parent ../Larvae_3_sorted.bam Chelymorpha_alternans.gff3 > count_table_Larvae_3_exon_Parent_new_ann.txt

