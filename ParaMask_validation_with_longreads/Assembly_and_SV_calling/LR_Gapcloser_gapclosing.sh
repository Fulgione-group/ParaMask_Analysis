#!/bin/bash
path_to_LR_Gapcloser="path/to/LR_Gapcloser/"
raw_data="/Path/to/raw_data.fastq"
assembly="/Path/to/Scaffolds.fasta"
working_dir="/Path/to/working_dir/"
output_dir="/Path/to/output_dir/"

cd ${working_dir}

#Index Scaffolds
samtools faidx ${assembly}

#Run LR Gapcloser
export PATH="${path_to_LR_Gapcloser}/src/:$PATH"
bash LR_Gapcloser.sh -i ${assembly} -l ${raw_data} -o ${output_dir}


