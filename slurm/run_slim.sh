#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

INHERITANCE=${1}
DOMINANCE=${2}
K=$3
BS=$4

slim -d g_size=999999 -d K=$K -d "r=log(2)" -d mig_rate=0.05 -d u=2.5e-8 -d rho=2.5e-8 -d b_s=$BS -d d_s=-0.0045 -d "inheritance='$INHERITANCE'" -d "dom_pattern='$DOMINANCE'" polyploid_all_DFE.slim 