#!/bin/bash
#$ -S /bin/bash
#$ -N modeltest
#$ -cwd


/ebio/ag-salem/projects/Metagenomic_assemblies/code/modeltest/bin/modeltest-ng -i yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed_alignment.phy -p 16 -d aa -T raxml


