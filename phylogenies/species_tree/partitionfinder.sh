#!/bin/bash
#$ -S /bin/bash
#$ -N pf_sym_nucl
#$ -cwd

set NUMEXPR_MAX_THREADS=128

source ~/.bashrc
conda activate py27
set NUMEXPR_MAX_THREADS=128

python PartitionFinder.py --raxml partitionfinder -p 128 --rcluster-max 100
