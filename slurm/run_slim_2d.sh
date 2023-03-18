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

### EXAMPLE: sbatch slurm/run_slim_2d.sh auto 10 recessive 50 0.000


INHERITANCE=${1}
D_SIZE=$2
DOMINANCE=${3}
K=$4
BS=$5

slim -d g_size=999999 -d K=$K -d d_size=$D_SIZE -d "r=log(2)" -d mig_rate=0.05 -d u=5e-8 -d rho=5e-8 -d b_s=$BS -d d_s=-0.001472 -d "inheritance='$INHERITANCE'" -d "dom_pattern='$DOMINANCE'" polyploid_all_DFE_2d.slim 