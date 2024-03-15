#!/bin/bash
for i in $(ls *mutation*_[02]_*SV*.dat);
do
	bash ~/scripts/ParaMask_simulation_scripts/transform_mutation_file.sh $i $(echo $i| rev| cut -d '.' -f2-| rev)_transformed.txt
done
