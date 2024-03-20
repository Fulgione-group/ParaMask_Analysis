#!/bin/bash

# Heterozygostiy/Genome Size Estimation - Jellyfish
ID="Accesion_ID"
path_raw_LR_data="/path/to/raw_data.fq"

jellyfish count -m 21 -s 100M -t 10 -C -F 2 <(zcat ${path_raw_LR_data}) -o ${ID}_mer_counts
#    -m: kmer length, 21 is commonly used
#    -s: size of hash table: should be genome size + extra kmers from seq errors
#    -t: number of threads
jellyfish histo -t 10 --high=1000000 ${ID}_mer_counts >${ID}_reads.histo			#Upload in GenomeScope
jellyfish dump mer_counts.jf > ${ID}_jf_dumps.fa

