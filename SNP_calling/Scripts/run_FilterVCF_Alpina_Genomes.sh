#!/bin/bash
java -jar ~/scripts/java/FilterVCF.jar --vcf /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATKcalls/Alpina_Genomes.vcf.gz\
	--indels-exclude\
	--biallelic\
	--quality 30\
	--depth 5,100000\
	--coverage /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATKcalls/GATK4.2_alpina_15_02_2023_chrall.vcf.avcov.txt,-1\
	--min-sample-cov 10\
	--gzip\
	--out /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATKcalls/Alpina_Genomes.filteredQ30LD5UD100K.vcf

