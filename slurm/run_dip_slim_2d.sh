#!/bin/bash
#SBATCH --partition=dschridelab
#SBATCH --constraint=rhel8
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

### EXAMPLE: sbatch slurm/run_dip_slim_2d.sh diff 0.5 10 recessive 50 0.000 0

DIP_MODEL=${1}
DIP_LAMBDA=${2}
D_SIZE=$3
DOMINANCE=${4}
K=$5
BS=$6
D_MUT=$7

slim -d g_size=999999 -d K=$K -d d_size=$D_SIZE -d "r=log(2)" -d dip_lambda=$DIP_LAMBDA -d remove_dip_muts=$D_MUT -d mig_rate=0.05 -d u_del=5e-8 -d u_dip=2.5e-10 -d u_ben=5e-9 -d rho=5e-8 -d b_s=$BS -d d_s=-0.001472 -d "dom_pattern='$DOMINANCE'" -d "dip_model='$DIP_MODEL'" polyploid_all_diploidization_2d.slim 
