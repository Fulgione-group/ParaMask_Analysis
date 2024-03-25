#!/bin/bash

ID="accession_ID"
raw_reads="/path/to/your/raw_reads.fastq"
reference="/path/to/your/reference.fa"        



# Map raw reads against reference
minimap2 -t 20 -ax map-hifi ${reference} ${raw_reads} --secondary=no | samtools sort -m 1G -o ${ID}_reads_against_reference.bam
