#!/bin/bash
#$ -S /bin/bash
#$ -N modeltest
#$ -cwd


modeltest-ng -i yellow_protein_alignment_outgroup.phy -p 16 -d aa -T raxml


