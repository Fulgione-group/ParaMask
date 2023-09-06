#!/bin/bash
awk 'BEGIN{FS=" ";OFS="\t"; printf"%s\t%s\t%s\t%s\n", "pos", "ac", "het" , "n"}{het=0;ac=0;for(i=2;i<=NF;i++){if($i==1){ac++}; if(i%2==0){a=$i}else{b=$i; if((a+b)==1){het++}; b=0;a=0}}; printf"%s\t%s\t%s\t%s\n", $1, ac, het , (NF-1)}' $1 > $2

