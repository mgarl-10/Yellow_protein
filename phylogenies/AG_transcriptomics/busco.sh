#!/bin/bash
#$ -S /bin/bash
#$ -N busco_Asparsa
#$ -cwd


source /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/etc/profile.d/conda.sh
conda activate busco

/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/bin/busco -m transcriptome -i Trinity.fasta -o busco_Acromis_sparsa_transcriptome -l /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/lineages/endopterygota_odb10 -f --offline --download_path /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/ -c 64

