#!/bin/bash

segments=$(ls mutations_new_1_Sim_SV_*_c0.05_paralogAll_tab.txt | wc -l)

head -n 1 mutations_new_1_Sim_SV_b1_l664_c0.05_paralogAll_tab.txt > Simulations_SV_SC_tab.txt

for (( i=1; i<=$segments; i++ ))
do

        p1=$(ls mutations_new_0_Sim_SC_b${i}_*_c0.05_tab.txt)
        p2=$(ls mutations_new_1_Sim_SV_b${i}_*_c0.05_paralogAll_tab.txt)
        tail -n +2 $p1 >> Simulations_SV_SC_tab.txt
        tail -n +2 $p2 >> Simulations_SV_SC_tab.txt
done

let segments+=1
p1=$(ls mutations_new_0_Sim_SC_b${segments}_*_c0.05_tab.txt)

tail -n +2 $p1 >> Simulations_SV_SC_tab.txt

