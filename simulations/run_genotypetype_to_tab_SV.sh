#!/bin/bash
for i in $(ls mutations_new_1_Sim_SV_*_c0.05_transformed.txt);
do
	bash genotypetype_to_tab_SV.sh $i $(echo $i| rev| cut -d '_' -f2-| rev)_paralog1_tab.txt
done

for i in $(ls mutations_new_2_Sim_SV_*_c0.05_transformed.txt);
do
        bash genotypetype_to_tab_SV.sh $i $(echo $i| rev| cut -d '_' -f2-| rev)_paralog2_tab.txt
done


