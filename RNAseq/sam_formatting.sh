#!/bin/bash
#$ -S /bin/bash
#$ -N sam
#$ -cwd

#samtools view -b -o Caplet_1.bam Caplet_1.sam
#samtools view -b -o Caplet_2.bam Caplet_2.sam
#samtools view -b -o Caplet_3.bam Caplet_3.sam
#samtools view -b -o Larvae_1.bam Larvae_1.sam
#samtools view -b -o Larvae_2.bam Larvae_2.sam
#samtools view -b -o Larvae_3.bam Larvae_3.sam
#samtools view -b -o Adult_FG_1.bam Adult_FG_1.sam
#samtools view -b -o Adult_FG_2.bam Adult_FG_2.sam
#samtools view -b -o Adult_FG_3.bam Adult_FG_3.sam
#samtools view -b -o Adult_OV_1.bam Adult_OV_1.sam
#samtools view -b -o Adult_OV_2.bam Adult_OV_2.sam
#samtools view -b -o Adult_OV_3.bam Adult_OV_3.sam


samtools sort -o Caplet_1_sorted.bam Caplet_1.bam
samtools sort -o Caplet_2_sorted.bam Caplet_2.bam
samtools sort -o Caplet_3_sorted.bam Caplet_3.bam
samtools sort -o Larvae_1_sorted.bam Larvae_1.bam
samtools sort -o Larvae_2_sorted.bam Larvae_2.bam
samtools sort -o Larvae_3_sorted.bam Larvae_3.bam
samtools sort -o Adult_FG_1_sorted.bam Adult_FG_1.bam
samtools sort -o Adult_FG_2_sorted.bam Adult_FG_2.bam
samtools sort -o Adult_FG_3_sorted.bam Adult_FG_3.bam
samtools sort -o Adult_OV_1_sorted.bam Adult_OV_1.bam
samtools sort -o Adult_OV_2_sorted.bam Adult_OV_2.bam
samtools sort -o Adult_OV_3_sorted.bam Adult_OV_3.bam
