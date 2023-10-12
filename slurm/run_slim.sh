#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o logs/slim.%A.out
#SBATCH -e logs/slim.%A.err

### EXAMPLE: sbatch slurm/run_slim.sh auto recessive 50 0.000 polyploid_all_DFE.slim

INHERITANCE=${1}
DOMINANCE=${2}
K=$3
BS=$4
RHO=$5
SCRIPT=$6

slim -d g_size=999999 -d K=$K -d "r=log(2)" -d mig_rate=0.005 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d b_s=$BS -d d_s=-0.005 -d rho=$RHO -d "inheritance='$INHERITANCE'" -d "dom_pattern='$DOMINANCE'" $SCRIPT  