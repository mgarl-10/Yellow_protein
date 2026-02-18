#!/bin/bash
#$ -S /bin/bash
#$ -N redun_C_alt
#$ -cwd


source ~/.bashrc
conda activate py27

python /code/redundans/redundans.py -i RNAseq_reads_1.fastq RNAseq_reads_2.fastq -l Chelymorpha_alternans_ccs.fasta -f Chelymorpha_alternans_hifiasm_assembly.fasta -o Chelymorpha_alternans_hifiasm_redundans -t 128 --log LOG
