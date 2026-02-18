#!/bin/bash
#$ -S /bin/bash
#$ -N blob
#$ -cwd

source ~/.bashrc
conda activate blobtools

blobtools create -i Chelymorpha_alternans_hifiasm_redun_assembly.fasta -b mapping_Chelymorpha_alternans.bam -t C_alternans_assembly_contaminants.txt -o blobplot



#/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/bin/blobtools/blobtools view -i blobplot.blobDB.json -o blobplot_C_bicolor

#/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/bin/blobtools/blobtools plot -i blobplot.blobDB.json -o plot 
/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/bin/blobtools/blobtools seqfilter -i Chelobasis_bicolor_hifiasm_assembly.fasta -l contigs_to_filter.txt -o C_bicolor_assembly_filtered -v 
