#!/bin/bash
#SBATCH --partition=dschridelab
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=16:00:00
#SBATCH -J plot
#SBATCH -o plot.%A.out
#SBATCH -e plot.%A.err

module load r
Rscript plot_revision.r