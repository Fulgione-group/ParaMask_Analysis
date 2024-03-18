#!/bin/bash

##### Fitting of coordinate systems ######
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5  ${ID}_assembly_only_chr.fasta ${reference} > ${ID}_reference_against_assembly.paf"
awk '{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9}' ${ID}_reference_against_assembly.paf > ${ID}_intermediate.txt
echo -e "query_name\tquery_start\tquery_end\treference\tref_start\tref_end" > Coordinates_${ID}_reference_purged.bed && cat ${ID}_intermediate.txt >> Coordinates_${ID}_reference_against_assembly.bed

#use Coordinates_${ID}_reference_against_assembly.bed as input for RStudio
