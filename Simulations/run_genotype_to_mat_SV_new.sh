#!/bin/bash
cd $2
for i in $(ls mutations_new_0_Sim_SV_*_c0.05_rep${1}_transformed.txt);
do
	par2=$(echo $i | cut -d'_' -f-2)"_2_"$(echo $i | cut -d'_' -f4-)
	bash ~/scripts/ParaMask_simulation_scripts/convert_dup_SV2.sh $i $par2 $(echo $i| rev| cut -d '_' -f2-| rev)_paralog3_mat.txt
	bash ~/scripts/ParaMask_simulation_scripts/genotype_to_mat_SV_new2.sh $(echo $i| rev| cut -d '_' -f2-| rev)_paralog3_mat.txt $i $(echo $i| rev| cut -d '_' -f2-| rev)_paralog1_mat.txt
        bash ~/scripts/ParaMask_simulation_scripts/genotype_to_mat_SV_new2.sh $(echo $i| rev| cut -d '_' -f2-| rev)_paralog3_mat.txt $par2 $(echo $i| rev| cut -d '_' -f2-| rev)_paralog2_mat.txt
done


