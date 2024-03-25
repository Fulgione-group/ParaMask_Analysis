#!/bin/bash

#Merge SV sets of sniffles and cuteSV

Survivor="/Path/to/survivor/"
ID="Accession_ID"

#Merge only overlapping SVs with allowed overlapping flanking region: 10% of SV size
${Survivor}/bin/SURVIVOR merge ${ID}_vcf_files.txt 0.1 2 1 0 1 35 ${ID}_only_overlapping_SVs.vcf


#Merge all SVs
${Survivor}/bin/SURVIVOR merge ${ID}_vcf_files.txt 0.1 1 1 0 1 35 ${ID}_all_SVs.vcf
