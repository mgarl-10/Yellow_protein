#!/bin/bash

for i in *yellow_species.fasta
do tag=${i%yellow_species.fasta}
/Users/mgarcialozano/code/muscle-5.1/src/Darwin/muscle -align "$i" -output "$tag"alignment.fasta
done 
