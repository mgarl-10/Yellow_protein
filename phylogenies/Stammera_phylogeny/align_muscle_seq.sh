#!/bin/bash

for i in *yellow_species.fasta
do tag=${i%yellow_species.fasta}
muscle -align "$i" -output "$tag"alignment.fasta
done 
