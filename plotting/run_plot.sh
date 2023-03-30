#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=4:00:00
#SBATCH -J plot
#SBATCH -o plot.%A.out
#SBATCH -e plot.%A.err

module load r
Rscript plot_file.r