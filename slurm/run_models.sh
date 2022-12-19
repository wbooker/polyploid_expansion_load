#run a series of models with different porameters
INHERITANCE=("diploid" "allo" "auto")
DOMINANCE=("additive" "overdominance_1" "overdominance_2" "overdominance_3")
K=(50)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for j in {1..5}
            do
            sbatch slurm/run_slim.sh $i $d $k 0.000
            done
        done
    done
done

INHERITANCE=("diploid" "auto" "allo")
DOMINANCE=("recessive" "underdominance_1" "underdominance_2")
K=(50)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for j in {1..5}
            do
            sbatch slurm/run_slim.sh $i $d $k 0.000
            done
        done
    done
done

INHERITANCE=("auto")
DOMINANCE=("duplex")
K=(50)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for j in {1..5}
            do
            sbatch slurm/run_slim.sh $i $d $k 0.000
            done
        done
    done
done

INHERITANCE=("diploid" "auto")
DOMINANCE=("br_dd" "bd_dr")
K=(50)
for i in "${INHERITANCE[@]}"
do
    for d in "${DOMINANCE[@]}"
    do
        for k in "${K[@]}"
        do
            for j in {1..5}
            do
            sbatch slurm/run_slim.sh $i $d $k 0.005
            done
        done
    done
done
