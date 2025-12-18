#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_HPC
#$ -cwd

raxml-ng-mpi -all -msa concatenated_aln_yellow_species.phy --prefix Stammera_yellow_species_protein --threads 16 -model partitions.txt --data-type AA --tree pars{25},rand{25} 


