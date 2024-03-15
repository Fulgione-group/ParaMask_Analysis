#!/bin/bash
##calculate according to Nei, hierfstat
awk 'BEGIN{FS=OFS="\t"; HoT=0; HeT=0}{if(NR>1){N++;Ho=$5; n=($2*$2);He=(n/(n-1))* (2*$4*(1-$4)-Ho/(2*n)); printf"%s\t%s\n", Ho, He; HoT+=Ho; HeT+=He}}END{printf"%s\t%s\t%s\t%s\n", "overall", HoT/N, HeT/N, 1-(HoT/HeT)}' GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES04Paramask.10PerMissing.snps.finalParaMaskIn.het.het.stat.txt > ES04_Fis.txt

