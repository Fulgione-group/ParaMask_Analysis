#!/bin/bash

counter=1
while IFS= read -r line
	do bash SeDuS_simSC_var.sh $line $counter 2 &
	counter=$(( counter +1 ))
        background=( $(jobs -p) )
        if (( ${#background[@]} == 40 )); then
                wait -n
        fi
done < SC_length_0.1_rep2.txt

