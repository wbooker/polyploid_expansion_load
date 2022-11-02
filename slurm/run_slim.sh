#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=96:00:00
#SBATCH --mem=16G
#SBATCH -J slim
#SBATCH -o slim.%A.out
#SBATCH -e slim.%A.err

slim -d K=50 -d "r=log(2)" -d mig_rate=0.05 -d u=2.5e-11 -d rho=2.5e-11 -d b_s=0.000 -d d_s=-0.005 -d "inheritance='auto'" -d "dom_pattern='recessive'" polyploid_all.slim 