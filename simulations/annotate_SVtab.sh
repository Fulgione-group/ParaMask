#!/bin/bash
awk 'BEGIN{FS=OFS="\t"; c=0;d=0}{if(NR==1){$4="block"}else{if((NR%2)==0){c++;$4="SC_b"c}else{d++;$4="SV_b"d}}; print $0}' SVtab.txt > SVtab_annotated.txt
