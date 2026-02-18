#!/bin/bash
#$ -S /bin/bash
#$ -N blob
#$ -cwd

source ~/.bashrc
conda activate blobtools

blobtools create -i Chelymorpha_alternans_hifiasm_redun_assembly.fasta -b mapping_Chelymorpha_alternans.bam -t C_alternans_assembly_contaminants.txt -o blobplot_Calternans

blobtools view -i blobplot.blobDB.json -o blobplot_Calternans

blobtools plot -i blobplot.blobDB.json -o blobplot_Calternans

blobtools seqfilter -i Chelymorpha_alternans_hifiasm_redun_assembly.fasta -l contigs_to_filter.txt -o C_alternans_assembly_filtered -v 
