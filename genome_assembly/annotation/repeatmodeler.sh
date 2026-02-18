#!/bin/bash
#$ -S /bin/bash
#$ -N rptm
#$ -cwd



#/ebio/ag-salem/projects/BeetleGenomes/code/RepeatModeler-2.0.3/BuildDatabase -name C_alternans -engine ncbi Chelymorpha_alternans_genome_filtered.fasta > database_build_C_alternans_run.out 

/ebio/ag-salem/projects/BeetleGenomes/code/RepeatModeler-2.0.3/RepeatModeler -database C_alternans -engine ncbi -pa 32 -LTRStruct > C_alternans_run_repeatmodeler.out 