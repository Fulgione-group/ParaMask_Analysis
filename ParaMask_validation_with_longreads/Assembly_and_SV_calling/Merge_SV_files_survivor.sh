#!/bin/bash

# SURVIVOR (merge both sets)
ID="ID"
/home/marimond/.conda/envs/bases/bin/SURVIVOR merge ${ID}_vcf_files.txt 0.1 2 1 0 1 35 ${ID}_only_overlapping_SVs.vcf
/home/marimond/.conda/envs/bases/bin/SURVIVOR merge ${ID}_vcf_files.txt 0.1 1 1 0 1 35 ${ID}_all_SVs.vcf
