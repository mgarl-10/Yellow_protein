#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_gene
#$ -cwd


raxml-ng-mpi -all -msa yellow_protein_alignment_outgroup.phy --prefix 7species_1outgroup_protein --threads 64 -model FLU+I+G4 --tree pars{25},rand{25} 


