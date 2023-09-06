#!/bin/bash
for i in $(ls *mutation*2_*SV*);
do
	bash transform_mutation_file.sh $i $(echo $i| rev| cut -d '.' -f2-| rev)_transformed.txt
done
