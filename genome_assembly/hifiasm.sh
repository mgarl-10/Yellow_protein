#!/bin/bash
#$ -S /bin/bash
#$ -N hifiasm
#$ -cwd


hifiasm -o Chelymorpha_alternans_hifiasm_assembly.asm -t 128 Chelymorpha_alternans_ccs.fasta
