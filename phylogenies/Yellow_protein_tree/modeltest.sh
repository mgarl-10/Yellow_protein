#!/bin/bash
#$ -S /bin/bash
#$ -N modeltest
#$ -cwd


modeltest-ng -i yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed_alignment.phy -p 16 -d aa -T raxml


