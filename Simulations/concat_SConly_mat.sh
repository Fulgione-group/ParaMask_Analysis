#!/bin/bash

cd $2
segments=$(ls mutations_new_1_Sim_SC_*_c0.05_rep${1}_mat.txt | wc -l)


for (( i=1; i<=$segments; i++ ))
do
        p1=$(ls mutations_new_1_Sim_SC_b${i}_*_c0.05_rep${1}_mat.txt)
        cat $p1 >> Simulations_SConly_rep${1}_mat.txt
done


