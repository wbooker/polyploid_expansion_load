#########run a series of models with different porameters

#INHERITANCE=("diploid" "allo" "auto")
#DOMINANCE=("DFE")
#K=(100)
#for i in "${INHERITANCE[@]}"
#do
#    for d in "${DOMINANCE[@]}"
#    do
#       for k in "${K[@]}"
#        do
#           for j in {1..25}
#            do
#            sbatch slurm/run_slim.sh $i $d $k 0.000 polyploid_all_DFE.slim
#           done
#       done
#   done
#done

#INHERITANCE=("diploid" "auto")
#DOMINANCE=("recessive")
#K=(100)
#b_s=(0.000 0.005)
#for i in "${INHERITANCE[@]}"
#do
#    for d in "${DOMINANCE[@]}"
#    do
#        for k in "${K[@]}"
#        do
#            for b in "${b_s[@]}"
##            do
#                for j in {1..25}
#                do
#                sbatch slurm/run_slim.sh $i $d $k $b polyploid_all_DFE.slim
#                done
#            done
#        done
#    done
#done

#INHERITANCE=("diploid" "auto")
#DOMINANCE=("bd_dr")
#K=(100)
#b_s=(0.005)
#for i in "${INHERITANCE[@]}"
#do
#    for d in "${DOMINANCE[@]}"
#    do
#        for k in "${K[@]}"
#        do
#            for b in "${b_s[@]}"
#            do
#                for j in {1..25}
#                do
#                sbatch slurm/run_slim.sh $i $d $k $b polyploid_all_DFE.slim
#                done
#            done
#        done
#    done
#done

INHERITANCE=("allo" "diploid")
DOMINANCE=("recessive" "additive" "DFE")
K=(100)
b_s=(0.000)
rho=(2.5e-8)
for i in "${INHERITANCE[@]}"
do
   for d in "${DOMINANCE[@]}"
   do
       for k in "${K[@]}"
       do
           for r in "${rho[@]}"
           do
               for j in {1..4}
               do
               sbatch slurm/run_slim.sh $i $d $k 0 $r polyploid_expansion.slim
               done
           done
       done
   done
done

# INHERITANCE=("diploid" "auto")
# DOMINANCE=("recessive" "additive")
# K=(100)
# b_s=(0.000)
# rho=(1e-7 1e-8 1e-9 1e-10 1e-6 1e-5)
# for i in "${INHERITANCE[@]}"
# do
#    for d in "${DOMINANCE[@]}"
#    do
#        for k in "${K[@]}"
#        do
#            for r in "${rho[@]}"
#            do
#                for j in {1..30}
#                do
#                sbatch slurm/run_slim.sh $i $d $k 0 $r polyploid_all_DFE.slim
#                done
#            done
#        done
#    done
# done

#INHERITANCE=("diploid" "auto" "allo")
#DOMINANCE=("additive")
#K=(100)
#b_s=(0.005)
#for i in "${INHERITANCE[@]}"
#do
###    for d in "${DOMINANCE[@]}"
#    do
#        for k in "${K[@]}"
#        do
#            for b in "${b_s[@]}"
#            do
#                for j in {1..25}
#                do
#                sbatch slurm/run_slim.sh $i $d $k $b polyploid_all_DFE.slim
#                done
#            done
#        done
#    done
#done

#DIP_LAMBDA=(100)
#DOMINANCE=("recessive")
#K=(100)
#d_s=(-0.01 -0.005 -0.0025)
#DIP_U=(1e-3 2e-3)
#D_MUT=(1)
#DEL_U=(2.5e-8)
#for i in "${DIP_LAMBDA[@]}"
#do
#   for d in "${DOMINANCE[@]}"
#    do
#         for k in "${K[@]}"
#         do
#             for s in "${d_s[@]}"
#             do
#                 for m in "${D_MUT[@]}"
#                 do
#                     for u in "${DIP_U[@]}"
#                     do
#                         for j in {1..28}
#                         do
#                         #sbatch slurm/run_dip_slim.sh diff $i $d $k 0.000 $s $u $m
#                         sbatch slurm/run_dip_slim.sh diff $i $d $k 0.000 $s $DEL_U $u 2.5e-8 $m 0 paper_run/pair_fit_dip_FIXED/test_recomb/no_peFit
#                         done
#                     done
#                 done
#             done
#         done
#     done
# done

#DIP_MODEL=(0.5 1.0 1.5 3.0)
#DOMINANCE=("recessive" "DFE")
#K=(100)
#d_s=(-0.005)
#D_MUT=(1 0)
#for i in "${DIP_MODEL[@]}"
#do
#    for d in "${DOMINANCE[@]}"
#    do
#        for k in "${K[@]}"
#        do
##            for s in "${d_s[@]}"
 #           do
 ##               for m in "${D_MUT[@]}"
##                do
#                    for j in {1..15}
#                    do
#                    sbatch slurm/run_dip_slim.sh diff $i $d $k 0.000 $s $m
#                    done
#                done
#            done
#        done
#    done
#done