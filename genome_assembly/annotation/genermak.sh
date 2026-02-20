#!/bin/bash
#$ -S /bin/bash
#$ -N gmark_Calt
#$ -cwd


source ~/.bashrc


/code/gmes_linux_64/gmes_petap.pl --ES --max_intron 460000 --soft_mask 2000 --cores 64 --sequence C_alternans_assembly.fa

