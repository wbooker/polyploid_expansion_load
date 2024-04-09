#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=72:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o logs/slim.%A.out
#SBATCH -e logs/slim.%A.err

### EXAMPLE: sbatch slurm/run_slim.sh auto recessive 50 -0.005 0.000 2.5e-8 polyploid_all_DFE.slim

INHERITANCE=${1}
DOMINANCE=${2}
K=$3
DS=$4
BS=$5
UD=$6
UB=$7
RHO=$8
SDIST=${9}
SCRIPT=${10}

slim -d "out_dir='output/revision_run/no_PAFC'" -d g_size=999999 -d K=$K -d "r=log(2)" -d mig_rate=0.05 -d u_del=$UD -d u_ben=$UB -d b_s=$BS -d d_s=$DS -d rho=$RHO -d "inheritance='$INHERITANCE'" -d "dom_pattern='$DOMINANCE'" -d "s_dist='$SDIST'" $SCRIPT  