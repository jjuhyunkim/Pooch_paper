#!/bin/bash

#SBATCH --partition=quick,norm
#SBATCH --cpus-per-task=10
#SBATCH --mem=50g
#SBATCH --array=0
#SBATCH --time=4:00:00
#SBATCH --output=logs/elasticnet.%j.log
#SBATCH --error=logs/elasticnet.%j.log
#SBATCH --job-name=elasticnet


export LD_LIBRARY_PATH=/data/kimj75/anaconda3/envs/plot_python/lib:$LD_LIBRARY_PATH

#chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
#chr=${chrs[$SLURM_ARRAY_TASK_ID]}

props=(5 10 15 20 25 30 40 50 60 70 80)
chr="chr22"
echo -e "$chr"
prop=${props[$SLURM_ARRAY_TASK_ID]}
echo -e "$prop"

python elasticnet.nestedcv.py $chr $prop
