#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J dip_slim
#SBATCH -o logs/dip_slim.%A.out
#SBATCH -e logs/dip_slim.%A.err

### EXAMPLE: sbatch slurm/run_dip_slim.sh diff 0.5 recessive 50 0.000 -0.005 0

DIP_MODEL=${1}
DIP_LAMBDA=${2}
DOMINANCE=${3}
K=$4
BS=$5
DS=$6
D_U=$7
D_MUT=$8 ##remove dip mutations prior to expansion (0 is don't remove)


slim -d g_size=999999 -d K=$K -d "r=log(2)" -d dip_lambda=$DIP_LAMBDA -d remove_dip_muts=$D_MUT -d mig_rate=0.05 -d u_del=2.5e-8 -d u_ben=2.5e-9 -d u_dip=$D_U -d b_s=$BS -d d_s=$DS -d "dom_pattern='$DOMINANCE'" -d "dip_model='$DIP_MODEL'" polyploid_all_diploidization_pairingFit_setLoci.slim 
