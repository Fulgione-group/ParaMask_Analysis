#!/bin/bash
VCF_file="/path/VCF_file.vcf"
ID="Sample_ID"

# Grep only duplications 
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\n' ${VCF_file} | grep 'DUP' | grep 'chr' > ${ID}_dup.txt
echo -e "Chromosome\tSTART\tSTOP\tLENGTH\tTYPE" > ${ID}_dup_with_header.txt && cat ${ID}_dup.txt >> ${ID}_dup_with_header.txt
