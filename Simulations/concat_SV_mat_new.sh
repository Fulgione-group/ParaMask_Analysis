#!/bin/bash
#how many
cd $2
segments=$(ls mutations_new_0_Sim_SV_*_c0.05_rep${1}_paralog2_mat.txt | wc -l)

for (( i=1; i<=$segments; i++ ))
do
	p1=$(ls mutations_new_0_Sim_SV_b${i}_*_c0.05_rep${1}_paralog1_mat.txt)
	p2=$(ls mutations_new_0_Sim_SV_b${i}_*_c0.05_rep${1}_paralog2_mat.txt)
	p3=$(ls mutations_new_0_Sim_SV_b${i}_*_c0.05_rep${1}_paralog3_mat.txt)
	out=$(echo $p1 | rev | cut -d'_' -f3-| rev)_paralogAll_mat.txt
	cat $p1 > tmp.txt
	cat $p2 >> tmp.txt
	if [[ $( head -1 $p3) != -1* ]]
	then
		cat $p3 >> tmp.txt
	fi
	sort -n -k1 tmp.txt >> $out
done
