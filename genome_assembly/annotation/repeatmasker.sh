#!/bin/bash
#$ -S /bin/bash
#$ -N rptm_Crubi
#$ -cwd



#/ebio/ag-salem/projects/BeetleGenomes/code/RepeatModeler-2.0.3/BuildDatabase -name C_rubiginosa -engine ncbi C_rubiginosa_genome_cleaned_sorted.fa > database_build_C_rubiginosa_run.out 

#/ebio/ag-salem/projects/BeetleGenomes/code/RepeatModeler-2.0.3/RepeatModeler -database C_rubiginosa -engine ncbi -pa 32 -LTRStruct > C_rubiginosa_run_repeatmodeler.out 

/ebio/ag-salem/projects/BeetleGenomes/code/maker/exe/RepeatMasker/RepeatMasker -pa 64 -gff -lib ./RM_106778.ThuFeb151425062024/consensi.fa.classified C_rubiginosa_genome_cleaned_sorted.fa > repeatmasker.out
