#! /usr/bin/env bash


./pal2nal.pl protein_aln_7species.fasta cds_7species.fasta -output fasta -nogap -codontable 1 > codon_based_aln_yellow_7species.fasta
