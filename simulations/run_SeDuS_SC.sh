#!/bin/bash

counter=1
while IFS= read -r line
	do bash run_SeDuS_simSC_var.sh $line $counter &
	counter=$(( counter +1 ))
        background=( $(jobs -p) )
        if (( ${#background[@]} == 20 )); then
                wait -n
        fi
done < SC_length_sample.txt

