#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J dip_slim
#SBATCH -o logs/dip_slim.%A.out
#SBATCH -e logs/dip_slim.%A.err

### EXAMPLE: sbatch slurm/run_dip_slim.sh diff 100 recessive 100 0.000 -0.005 2.5e-8 2.5e-9 1e-3 1e-6 fixed 1 1 85 1 output/test_dip
DIP_MODEL=${1}
DIP_LAMBDA=${2}
DOMINANCE=${3}
K=$4
BS=$5
DS=$6
DEL_U=$7
BEN_U=$8
DIP_U=$9
RHO=${10}
SDIST=${11}
D_MUT=${12} ##remove dip mutations prior to expansion (0 is don't remove)
OUT_DIR=${13}

slim -d "out_dir='$OUT_DIR'" -d g_size=999999 -d K=$K -d "r=log(2)" -d dip_lambda=$DIP_LAMBDA -d rho=$RHO -d remove_dip_muts=$D_MUT -d mig_rate=0.05 -d u_del=$DEL_U -d u_ben=$BEN_U -d u_dip=$DIP_U -d b_s=$BS -d d_s=$DS -d "dom_pattern='$DOMINANCE'" -d "dip_model='$DIP_MODEL'" -d "s_dist='$SDIST'" polyploid_diploidization_dominant_meiotic.slim