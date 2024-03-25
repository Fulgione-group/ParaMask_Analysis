#!/bin/bash
java -jar /home/btjeng/scripts/java/PrepareParaMaskInput_fromVCF.jar \
-VF \
--vcf /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATKfinal/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.b.vcf.gz \
--out /netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATKfinal/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.finalParaMaskIn.het \
--missingness 0.1

