#!/bin/bash
for i in $(ls mutations_new_0_Sim_SC_*_c0.05_transformed.txt);
do
	bash genotypetype_to_tab_SC.sh $i $(echo $i| rev| cut -d '_' -f2-| rev)_tab.txt
done
