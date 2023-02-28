#!/bin/bash 
touch core_fits_2501.txt
for i in {1..50}
do
    cat $i/${i}_log.csv | head -n 2 | tail -n 1 | awk -F, '{print $2}' >> core_fits_2501.txt
done