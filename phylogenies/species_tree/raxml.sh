#!/bin/bash
#$ -S /bin/bash
#$ -N raxml_species
#$ -cwd


raxml-ng-mpi -all -msa concatenated_alignment_7species.phy --prefix host_7species_unrooted --threads 16 -model partitions.txt --tree pars{25},rand{25} 


