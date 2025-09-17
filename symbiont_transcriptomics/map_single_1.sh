#!/bin/bash
#$ -S /bin/bash
#$ -N map1_sing
#$ -cwd

source ~/.bashrc
conda activate bowtie2

bowtie2 -x Stammera -U CTL_UNT_R1_trimmed_single_1.fastq.gz --fr --no-unal -p 64 -S CTL_UNT_R1_Stammera_single.sam

bowtie2 -x Stammera -U CTL_HUM_R1_trimmed_single_1.fastq.gz --fr --no-unal -p 64 -S CTL_HUM_R1_Stammera_single.sam