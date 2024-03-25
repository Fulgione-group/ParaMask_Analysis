#!/bin/bash

ID="Accession_ID"
reference="/Path/to/reference.fasta"

# Map Assembly against reference
minimap2 -t 20 -cx asm5  ${ID}_assembly_only_chr.fasta ${reference} > ${ID}_reference_against_assembly.paf

# Process PAF file (get coordinates of query and target)
awk '{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $8 "\t" $9}' ${ID}_reference_against_assembly.paf > ${ID}_intermediate.txt
echo -e "query_name\tquery_start\tquery_end\treference\tref_start\tref_end" > Coordinates_${ID}_reference_purged.bed && cat ${ID}_intermediate.txt >> Coordinates_${ID}_reference_against_assembly.bed

#use Coordinates_${ID}_reference_against_assembly.bed as input for RStudio
