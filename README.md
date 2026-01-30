

Bioinformatics analysis pipeline for: 

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H. Yellow protein co-opted to sustain obligate symbiosis in beetles.

# 1. Comparative transcriptomics 

cd comparative_transcriptomics

## 1.1 Quality control

Requires raw fastq.gz paired-end reads in this directory as well as TruSeq3-PE.fa file with adapter sequences

1. Remove adapters and quality trimming: qsub trimmomatic.sh

## 1.2 Map quality trimmed reads to reference genome 

1. Index reference genome: 
2. Map trimmed reads to _Chelymorpha alternans_ genome by Hisat2:
3. Convert bam files to sorted sam files with samtools: 

## 1.3 Quantify mapped reads



