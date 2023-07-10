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

#INHERITANCE=("auto")
#DOMINANCE=("duplex")
#K=(100)
#b_s=(0.000)
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

DIP_LAMBDA=(100)
DOMINANCE=("recessive")
K=(100)
d_s=(-0.005 0.01 0.0025)
DIP_U=(0.001 0.002)
D_MUT=(0)
for i in "${DIP_LAMBDA[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for s in "${d_s[@]}"
            do
                for m in "${D_MUT[@]}"
                do
                    for u in "${DIP_U[@]}"
                    do
                        for j in {1..25}
                        do
                        sbatch slurm/run_dip_slim.sh diff $i $d $k 0.000 $s $u $m
                        done
                    done
                done
            done
        done
    done
done

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