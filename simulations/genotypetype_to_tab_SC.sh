#!/bin/bash

##Error rate for heterozygotes is 0.01 either homozygous
awk 'BEGIN{
	FS=" ";
	OFS="\t";
	printf"%s\t%s\t%s\t%s\n", "pos", "ac", "het" , "n"
}{
	het=0;
	ac=0;
	for(i=2;i<=NF;i++){
		if(i%2==0){
			a=$i
		}else{
			b=$i;
			#homozygous derived 0.98 prob to appear heterozygous, 0.01 prob to appear homozygous ancestral or derived 
			if((a+b)==1){
				r=rand();
				if(r<0.01){
					ac=(ac+2)
				}else if(r>0.02){
					het++;
					ac++;
				}
			}
			if((a+b)==2){
				ac=(ac +2);
			}
			b=0;
			a=0
		}
	};
	printf"%s\t%s\t%s\t%s\n", $1, ac, het , (NF-1)
}' $1 > $2

