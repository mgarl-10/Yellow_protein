#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_species
#$ -cwd


/ebio/ag-salem/projects/Metagenomic_assemblies/code/raxml-ng_v1.2.0_linux_x86_64_MPI/bin/raxml-ng-mpi -all -msa concatenated_alignment_7species.phy --prefix host_7species_unrooted --threads 16 -model partitions.txt --tree pars{25},rand{25} 


