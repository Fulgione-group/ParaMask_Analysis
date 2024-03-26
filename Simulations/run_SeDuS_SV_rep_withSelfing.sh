#!/bin/bash

counter=1
while IFS= read -r line
	do bash SeDuS_simSV_var_withSelfing.sh $line $counter REP_NUMBER &
	counter=$(( counter +1 ))
        background=( $(jobs -p) )
        if (( ${#background[@]} == 60 )); then
                wait -n
        fi
done < SV_length_0.5_rep1.txt

