#!/bin/bash

table_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/assemblies_for_script.txt"
while IFS=$'\t' read -r col1 col2
do

ID="${col1}"
raw_data="${col2}"
genome_size="ESTIMATED_GENOME_SIZE"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye/
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi ${raw_data} -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye -g ${genome_size} -t 45"
done < "$table_file"
