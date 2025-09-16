#!/bin/bash
#$ -S /bin/bash
#$ -N pf_sym_prot
#$ -cwd

set NUMEXPR_MAX_THREADS=128

source ~/.bashrc
conda activate py27
set NUMEXPR_MAX_THREADS=128

python /ebio/ag-salem/projects/Metagenomic_assemblies/code/partitionfinder-2.1.1/PartitionFinderProtein.py --raxml partitionfinder -p 128 --rcluster-max 100
