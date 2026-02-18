#!/bin/bash
#$ -S /bin/bash
#$ -N fun_Crubi
#$ -cwd


source ~/.bashrc

conda activate funannotate_2

#funannotate clean -i C_rubiginosa_assembly_filtered.fasta --minlen 1000 -o C_rubiginosa_genome_cleaned.fa --cpus 32
#funannotate sort -i C_rubiginosa_genome_cleaned.fa -b contig -o C_rubiginosa_genome_cleaned_sorted.fa 

#funannotate mask -i C_rubiginosa_genome_cleaned_sorted.fa -m repeatmasker -l ./RM_95950.WedFeb211359442024/consensi.fa.classified --cpus 64 -o C_rubiginosa_assembly.fa

#funannotate train -i C_rubiginosa_assembly.fa -o fun --left /ebio/ag-salem/projects/CassidinaeGenomics/Marleny_Crubiginosa/RNA_S15195Nr9.1.fastq.gz --right /ebio/ag-salem/projects/CassidinaeGenomics/Marleny_Crubiginosa/RNA_S15195Nr9.2.fastq.gz --stranded RF --species "Cassida rubiginosa" --cpus 32 --max_intronlen 460000 --memory 100G


#funannotate predict -i C_rubiginosa_assembly.fa -o fun -s "Cassida rubiginosa" --max_intronlen 460000 --organism other --repeats2evm --busco_db endopterygota --weights glimmerhmm:0 --genemark_gtf /ebio/ag-salem/projects/CassidinaeGenomics/data/C_rubiginosa_Annotation/funannotate/genemark/genemark.gtf


#funannotate update -i fun  --cpus 64

#funannotate species -s cassida_rubiginosa -a fun/predict_results/cassida_rubiginosa.parameters.json
#funannotate iprscan -i fun -m local --iprscan_path /ebio/ag-salem/projects/BeetleGenomes_Annotation/code/my_interproscan/interproscan-5.55-88.0/interproscan.sh -c 64
funannotate annotate -i fun --species "Cassida rubiginosa" --busco_db endopterygota --cpus 64
