#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate baconda

vcftools --gzvcf GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03ES04Paramask.10PerMissing.snps.0.1.b.vcf.gz\
 --out GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03ES04Paramask.10PerMissing.snps.0.1.b.vcf.FST\
 --weir-fst-pop ES03_1000Genomes_IDonly.txt\
 --weir-fst-pop ES04_1000Genomes_IDonly.txt\
 --fst-window-size 400000000\  ##big window to cover the genome
 --fst-window-step 400000000
wait
conda deactivate
