#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_gene
#$ -cwd


/ebio/ag-salem/projects/Metagenomic_assemblies/code/raxml-ng_v1.2.0_linux_x86_64_MPI/bin/raxml-ng-mpi -all -msa yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed_alignment.phy --prefix yellow_all --threads 16 -model LG+G4


