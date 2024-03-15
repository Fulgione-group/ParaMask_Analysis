#!/bin/bash

cd $2
segments=$(ls mutations_new_0_Sim_SV_*_c0.05_rep${1}_paralogAll_mat.txt | wc -l)


for (( i=1; i<=$segments; i++ ))
do
        p1=$(ls mutations_new_1_Sim_SC_b${i}_*_c0.05_rep${1}_mat.txt)
        p2=$(ls mutations_new_0_Sim_SV_b${i}_*_c0.05_rep${1}_paralogAll_mat.txt)
        cat $p1 >> Simulations_SV_SC_rep${1}_mat.txt
        cat $p2 >> Simulations_SV_SC_rep${1}_mat.txt
done

let segments+=1
p1=$(ls mutations_new_1_Sim_SC_b${segments}_*_c0.05_rep${1}_mat.txt)

cat $p1 >> Simulations_SV_SC_rep${1}_mat.txt

