#!/bin/bash
java -jar /home/btjeng/scripts/java/PrepareParaMaskInput_fromVCF.jar \
-VF \
--vcf /GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.b.vcf.gz \
--out /GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.finalParaMaskIn.het \
--missingness 0.1

