#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8 
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

slim pieschl_sim.slim