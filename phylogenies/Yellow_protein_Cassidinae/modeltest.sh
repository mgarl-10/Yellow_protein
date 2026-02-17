#!/bin/bash
#$ -S /bin/bash
#$ -N modeltest
#$ -cwd


/ebio/ag-salem/projects/Metagenomic_assemblies/code/modeltest/bin/modeltest-ng -i yellow_protein_alignment_outgroup.phy -p 16 -d aa -T raxml


