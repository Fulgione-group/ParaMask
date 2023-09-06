#!/bin/bash
#how many
segments=$(ls mutations_new_2_Sim_SV_*_c0.05_paralog2_tab.txt | wc -l)

for (( i=1; i<=$segments; i++ ))
do
	p1=$(ls mutations_new_1_Sim_SV_b${i}_*_c0.05_paralog1_tab.txt)
	p2=$(ls mutations_new_2_Sim_SV_b${i}_*_c0.05_paralog2_tab.txt)
	out=$(echo $p1 | rev | cut -d'_' -f3-| rev)_paralogAll_tab.txt
	tail -n +2 $p1 > tmp.txt
	tail -n +2 $p2 >> tmp.txt
	head -n 1 $p1 > $out
	sort -n -k1 tmp.txt >> $out

done
