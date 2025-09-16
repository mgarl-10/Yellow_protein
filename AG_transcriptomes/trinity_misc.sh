#!/bin/bash
#$ -S /bin/bash
#$ -N trinity_abun_Acromis
#$ -cwd

perl /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/funannotate_2/bin/TrinityStats.pl Trinity.fasta > trinityStats.log
perl /ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/funannotate_2/bin/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left ../Acromis_sparsa_AG_trimmed_1.fastq.gz --right ../Acromis_sparsa_AG_trimmed_2.fastq.gz --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir_Acromis_sparsa --coordsort_bam
