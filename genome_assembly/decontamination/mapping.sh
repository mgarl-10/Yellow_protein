#!/bin/bash
#$ -S /bin/bash
#$ -N bwa
#$ -cwd

bwa index Chelymorpha_alternans_hifiasm_redun_assembly.fasta

bwa mem -t 64 Chelymorpha_alternans_hifiasm_redun_assembly.fasta Chelymorpha_alternans_ccs.fasta > mapping_Chelymorpha_alternans.sam

samtools sort --threads 64 mapping_Chelymorpha_alternans.sam > mapping_Chelymorpha_alternans.bam 

samtools index mapping_Chelymorpha_alternans.bam 
