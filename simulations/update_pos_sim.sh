#!/bin/bash
awk 'BEGIN{FS=OFS="\t"}{if(NR==FNR){if(FNR>1){p[$4]=($1-1)}}else{if(FNR>1){$1=($1+p[$5])}; print $0}}' SVtab_annotated.txt Simulations_SV_SC_tab.txt > Simulations_SV_SC_tab_cpos.txt
