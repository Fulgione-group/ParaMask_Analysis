#!/bin/bash

counter=1
while IFS= read -r line
	do bash SeDuS_simSV_var.sh $line $counter 1 &
	counter=$(( counter +1 ))
        background=( $(jobs -p) )
        if (( ${#background[@]} == 40 )); then
                wait -n
        fi
done < SV_length_0.1_rep1.txt

