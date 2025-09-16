#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_HPC
#$ -cwd


/ebio/ag-salem/projects/Metagenomic_assemblies/code/raxml-ng_v1.2.0_linux_x86_64_MPI/bin/raxml-ng-mpi -all -msa concatenated_aln_yellow_species.phy --prefix Stammera_yellow_species_protein --threads 16 -model partitions.txt --data-type AA --tree pars{25},rand{25} 


