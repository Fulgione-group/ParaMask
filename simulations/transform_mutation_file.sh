#!/bin/bash
awk 'BEGIN{FS=" "; OFS="\t"}{if(NR==1){ split(FILENAME, c , "_"); split(c[7], d , "l");l=d[2]}else if(NR==5){for(i=2; i<=(NF-1);i++){printf"%s\t", int($i*l)}; printf"%s\n", int($NF*l)}else if(NR>5){split($0,a, "");for(i=1; i<=(length(a)-1);i++){printf"%s\t", a[i]}; printf"%s\n", a[length(a)]}}' $1 | awk 'BEGIN{FS=OFS="\t"}{for (i=1; i<=NF; i++){a[NR,i] = $i}}NF>p{ p = NF }END {for(j=1; j<=p; j++){str=a[1,j];for(i=2; i<=NR; i++){str=str" "a[i,j]}print str}}' > $2


