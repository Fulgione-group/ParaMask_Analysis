#!/bin/bash

counter=1
while IFS= read -r line
	do bash SeDuS_simSC_var_withSelfing.sh $line $counter REP_NUMBER &
	counter=$(( counter +1 ))
        background=( $(jobs -p) )
        if (( ${#background[@]} == 60 )); then
                wait -n
        fi
done < SC_length_0.5_rep1.txt

