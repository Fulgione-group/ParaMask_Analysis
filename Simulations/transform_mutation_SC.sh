#!/bin/bash
cd $1
for i in $(ls *mutation*1_*SC*.dat);
do
	bash ~/scripts/ParaMask_simulation_scripts/transform_mutation_file.sh $i $(echo $i| rev| cut -d '.' -f2-| rev)_transformed.txt
done
