

Bioinformatics analysis pipeline for: 

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H. Yellow protein co-opted to sustain obligate symbiosis in beetles.

# 1. Comparative transcriptomics 

cd comparative_transcriptomics

## 1.1 Quality control

cd comparative_transcriptomics/quality_control

Requires raw fastq.gz paired-end reads in this directory as well as TruSeq3-PE.fa file with adapter sequences

1. Remove adapters and quality trimming: qsub trimmomatic.sh
2. Check read quality: qsub fastqc.sh 

## 1.2 Map quality trimmed reads to reference genome 

cd comparative_transcriptomics/mapping

1. Index reference genome: 
2. Map trimmed reads to _Chelymorpha alternans_ genome:
3. Convert bam files to sorted sam files:
4. 

