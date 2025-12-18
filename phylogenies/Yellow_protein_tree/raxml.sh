#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_gene
#$ -cwd


raxml-ng-mpi -all -msa yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed_alignment.phy --prefix yellow_all --threads 16 -model LG+G4


