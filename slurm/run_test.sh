#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

slim test_dip_fitEffect.slim 