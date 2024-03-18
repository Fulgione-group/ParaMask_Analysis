#!/bin/bash

work_dir="/path/to/your/workingdirectory"

ID="accession_ID"

raw_reads="/path/to/your/raw_reads.fastq"

reference="/path/to/your/reference.fa"        #/netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta

# set working directory
cd ${work_dir}


# Map raw reads against reference
minimap2 -t 20 -ax map-hifi ${reference} ${raw_reads} --secondary=no | samtools sort -m 1G -o ${ID}_reads_against_reference.bam
