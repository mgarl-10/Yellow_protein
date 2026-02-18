#!/bin/bash
#$ -S /bin/bash
#$ -N fun_Calt
#$ -cwd


source ~/.bashrc

conda activate funannotate


funannotate train -i C_alternans_assembly.fa -o fun --left RNAseq_reads.1.fastq.gz --right RNAseq_reads.1.fastq.gz --stranded RF --species "Chelymorpha alternans" --cpus 32 --max_intronlen 460000 --memory 100G


funannotate predict -i C_alternans_assembly.fa -o fun -s "Chelymorpha_alternans" --max_intronlen 460000 --organism other --repeats2evm --busco_db endopterygota --weights glimmerhmm:0 --genemark_gtf genemark.gtf


funannotate update -i fun  --cpus 64

funannotate species -s chelymorpha_alternans -a fun/predict_results/cassida_rubiginosa.parameters.json

funannotate iprscan -i fun -m local --iprscan_path /code/interproscan/interproscan.sh -c 64

funannotate annotate -i fun --species "Chelymorpha alternans" --busco_db endopterygota --cpus 64
