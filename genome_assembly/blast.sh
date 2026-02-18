#!/bin/bash
#$ -S /bin/bash
#$ -N blob
#$ -cwd


blastn -query Chelymorpha_alternans_hifiasm_assembly.fasta -db nt -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 1 -max_hsps 1 -num_threads 64 -evalue 1e-25 -out C_alternans_assembly_contaminants.txt
