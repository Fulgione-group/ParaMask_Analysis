#!/bin/bash

table_file="/Path/to/table_file.txt"  #including accession ID (col1) and raw read paths (col 2)
  while IFS=$'\t' read -r col1 col2
    do

    ID="${col1}"
    raw_data="${col2}"
    genome_size="ESTIMATED_GENOME_SIZE"
    working_dir="/Path/tp/work_dir/"
    output_dir="/Path/to/output_dir/"

    cd ${working_dir}
    bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi ${raw_data} -o ${output_dir} -g ${genome_size} -t 45"
done < "$table_file"
