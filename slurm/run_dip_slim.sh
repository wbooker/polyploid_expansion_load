#!/bin/bash
#SBATCH --partition=general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH -J dip_slim
#SBATCH -o logs/dip_slim.%A.out
#SBATCH -e logs/dip_slim.%A.err

### EXAMPLE: sbatch slurm/run_dip_slim.sh diff 100 recessive 100 0.000 -0.000 2.5e-10 1e-3 2.5e-8 1 1 paper_run/pair_fit_dip_FIXED/test_recomb/no_s

DIP_MODEL=${1}
DIP_LAMBDA=${2}
DOMINANCE=${3}
K=$4
BS=$5
DS=$6
DEL_U=$7
DIP_U=$8
RHO=$9
D_MUT=${10} ##remove dip mutations prior to expansion (0 is don't remove)
M_F=${11} ## Fitness cost to pairing efficiency? (1 is yes, 0 no)
OUT_DIR=${12}


slim -d "out_dir='$OUT_DIR'" -d g_size=999999 -d K=$K -d "r=log(2)" -d dip_lambda=$DIP_LAMBDA -d rho=$RHO -d remove_dip_muts=$D_MUT -d meiotic_fitness=$M_F -d mig_rate=0.05 -d u_del=$DEL_U -d u_ben=2.5e-9 -d u_dip=$DIP_U -d b_s=$BS -d d_s=$DS -d "dom_pattern='$DOMINANCE'" -d "dip_model='$DIP_MODEL'" polyploid_all_diploidization_pairingFit_setLoci.slim