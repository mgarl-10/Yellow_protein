#!/bin/bash
#$ -S /bin/bash
#$ -N rptm_Calt
#$ -cwd



/code/RepeatModeler-2.0.3/BuildDatabase -name C_alternans -engine ncbi C_alternans_genome_cleaned_sorted.fa > database_build_C_alternans_run.out 

/code/RepeatModeler-2.0.3/RepeatModeler -database C_alternans -engine ncbi -pa 32 -LTRStruct > C_alternans_run_repeatmodeler.out 

/code/RepeatMasker/RepeatMasker -pa 64 -gff -lib ./RM/consensi.fa.classified C_alternans_genome_cleaned_sorted.fa > repeatmasker.out
