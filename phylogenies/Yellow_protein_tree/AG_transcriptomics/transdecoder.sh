#!/bin/bash
#$ -S /bin/bash
#$ -N transdec_Acromis
#$ -cwd


#/ebio/ag-salem/projects/BeetleGenomes_Annotation/code/TransDecoder/TransDecoder.LongOrfs -t Acromis_sparsa_AG_trinity_out/Trinity.fasta --output_dir Acromis_sparsa_transdecoder

#hmmscan --cpu 16 --domtblout pfam_Acromis_sparsa.domtblout /ebio/ag-salem/projects/BeetleGenomes_Annotation/code/Trinotate/Pfam-A.hmm ./Acromis_sparsa_transdecoder/longest_orfs.pep


/ebio/ag-salem/projects/BeetleGenomes_Annotation/code/TransDecoder/TransDecoder.Predict -t Acromis_sparsa_AG_trinity_out/Trinity.fasta --retain_pfam_hits pfam_Acromis_sparsa.domtblout --cpu 16 --output_dir Acromis_sparsa_transdecoder
