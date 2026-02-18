#!/bin/bash
#$ -S /bin/bash
#$ -N busco_hifi
#$ -cwd


source /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/etc/profile.d/conda.sh
conda activate busco

/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/bin/busco -m genome -i contigs.reduced.fa -o busco_redundans_hifiasm -l /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/lineages/endopterygota_odb10 -f --offline --download_path /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/ -c 64

