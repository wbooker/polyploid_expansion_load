#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

slim -seed 16948153676115 polyploid_allo_noMatch.slim