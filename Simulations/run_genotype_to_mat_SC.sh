#!/bin/bash
cd $2
for i in $(ls mutations_new_1_Sim_SC_*_c0.05_rep${1}_transformed.txt);
do
	bash ~/scripts/ParaMask_simulation_scripts/genotype_to_mat_SC_new.sh $i $(echo $i| rev| cut -d '_' -f2-| rev)_mat.txt
done
