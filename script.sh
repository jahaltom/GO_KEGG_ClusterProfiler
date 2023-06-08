#!/bin/bash
#SBATCH -N 1
#SBATCH -p defq
#SBATCH --ntasks-per-node 24
#SBATCH -t 00:20:00


module load R

Rscript script.r
