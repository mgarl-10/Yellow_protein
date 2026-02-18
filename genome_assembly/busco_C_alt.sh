#!/bin/bash
#$ -S /bin/bash
#$ -N busco_hifi
#$ -cwd


source ~/.bashrc
conda activate busco

busco -m genome -i Chelymorpha_alternans_hifiasm_redun.fasta -o busco_redundans_hifiasm -l endopterygota_odb10 -f --offline --download_path /busco/busco_downloads/ -c 64

