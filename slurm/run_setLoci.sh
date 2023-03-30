#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

### EXAMPLE: sbatch slurm/run_dip_slim.sh diff 0.5 recessive 50 0.000 0

DIP_MODEL=${1}
DIP_LAMBDA=${2}
DOMINANCE=${3}
K=$4
BS=$5
D_MUT=$6 ##remove dip mutations prior to expansion (0 is don't remove)


slim -d g_size=999999 -d K=$K -d "r=log(2)" -d dip_lambda=$DIP_LAMBDA -d remove_dip_muts=$D_MUT -d mig_rate=0.05 -d u_r=5e-8 -d u_d=1e-4 -d u_b=5e-8 -d rho=5e-8 -d b_s=$BS -d d_s=-0.001472 -d "dom_pattern='$DOMINANCE'" -d "dip_model='$DIP_MODEL'" polyploid_all_diploidization_setLoci.slim 
