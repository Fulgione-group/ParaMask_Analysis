#!/bin/bash

java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03Paramask.10PerMissing.snps.b.vcf.gz\
        --cutoff 0.1\
        --seqlength 171043752\ #accounting for density with missing data of all site VCF
        --id ES03_collapsed\
        --replicate 1

