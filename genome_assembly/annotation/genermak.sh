#!/bin/bash
#$ -S /bin/bash
#$ -N gmark_Crubi
#$ -cwd


source ~/.bashrc


/ebio/ag-salem/projects/BeetleGenomes_Annotation/code/maker/exe/gmes_linux_64/gmes_petap.pl --ES --max_intron 460000 --soft_mask 2000 --cores 64 --sequence ../C_rubiginosa_assembly.fa

