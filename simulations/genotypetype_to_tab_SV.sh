#!/bin/bash

##Error rate for heterozygotes is 0.01 either homozygous
awk 'BEGIN{
	FS=" ";
	OFS="\t";
	printf"%s\t%s\t%s\t%s\t%s\n", "pos", "ac", "het" , "n", "block"
}{
	het=0;
	ac=0;
	n=0;
	split(FILENAME, z,"_");
	block=z[6]
	for(i=2;i<=NF;i++){
		if(i%2==0){
			a=$i
		}else{
			b=$i;
			#if heterozygous 0.2 prob to be collapsed heterozygous, 0.2 prob to be homozygous ancestral and 0.6 prob to be missing 
			if((a+b)==1){
				r=rand();
				if(r<0.2){
					het++;
					ac++;
					n=(n+2)
				}else if(r<0.4){
					n=(n+2)
				}
			#homozygous derived 0.98 prob to appear heterozygous, 0.01 prob to appear homozygous ancestral or derived 
			}else if((a+b)==2){
				r=rand();
                                if(r<0.01){
                                        ac=(ac+2)
					n=(n+2)
				}else if(r< 0.02){
					n=(n+2)
                                }else{
                                        het++;
                                        ac++;
					n=(n+2)
                                }
			}else{
				n=(n+2)
			}
			b=0;
			a=0;
		}
	};
	if(ac>0){
		printf"%s\t%s\t%s\t%s\t%s\n", $1, ac, het , n, "SV_"block;
	}
}' $1 > $2
