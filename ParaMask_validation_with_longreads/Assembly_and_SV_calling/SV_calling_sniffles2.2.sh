#!/bin/bash
# SV calling using sniffles v2.2 (as conda env)

conda activate sniffles_env

input_bam="/path/to/input/bam" 
ID="accession_ID"
repeat_file="/Path/to/repeats.bed"
reference="/Path/to/reference.fasta"

sniffles --input input_bam --snf ${ID}_sniffles2_2_repeats.snf --threads 20 --vcf ${ID}_sniffles2.2_repeats.vcf --reference ${reference} --threads 20 --minsvlen 35 --tandem-repeats ${repeat_file} --allow-overwrite
