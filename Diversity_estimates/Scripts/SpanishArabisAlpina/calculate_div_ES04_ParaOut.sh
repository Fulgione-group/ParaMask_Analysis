#!/bin/bash

java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.b.vcf_ParaOut.vcf.gz\
        --cutoff 0.1\
        --seqlength 266128384\
        --id ES04_ParaOut\
        --replicate 1

