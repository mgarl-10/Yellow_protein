#!/bin/bash


for i in *_renamed.fa
do tag=${i%*_renamed.fa}
./faSomeRecords.py -f "$i" --list species.txt --out "$tag"_yellow_species.fasta
done
