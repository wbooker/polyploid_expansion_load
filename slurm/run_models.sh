#run a series of models with different porameters
INHERITANCE=("diploid" "allo" "auto")
DOMINANCE=("DFE")
K=(50)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for j in {1..50}
            do
            sbatch slurm/run_slim.sh $i $d $k 0.000
            done
        done
    done
done

INHERITANCE=("diploid" "auto")
DOMINANCE=("recessive")
K=(50)
b_s=(0.000 0.0045)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for b in "${b_s[@]}"
            do
                for j in {1..50}
                do
                sbatch slurm/run_slim.sh $i $d $k $b
                done
            done
        done
    done
done

INHERITANCE=("diploid" "auto" "allo")
DOMINANCE=("additive")
K=(50)
b_s=(0.000 0.0045)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for b in "${b_s[@]}"
            do
                for j in {1..50}
                do
                sbatch slurm/run_slim.sh $i $d $k $b
                done
            done
        done
    done
done
