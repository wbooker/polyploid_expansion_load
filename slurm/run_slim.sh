#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

### EXAMPLE: sbatch slurm/run_slim.sh auto recessive 50 0.000

INHERITANCE=${1}
DOMINANCE=${2}
K=$3
BS=$4

slim -d g_size=999999999 -d K=$K -d "r=log(2)" -d mig_rate=0.05 -d u_del=2.5e-11 -d u_ben=2.5e-11 -d rho=2.5e-11 -d b_s=$BS -d d_s=-0.001472 -d "inheritance='$INHERITANCE'" -d "dom_pattern='$DOMINANCE'" polyploid_all_DFE.slim 