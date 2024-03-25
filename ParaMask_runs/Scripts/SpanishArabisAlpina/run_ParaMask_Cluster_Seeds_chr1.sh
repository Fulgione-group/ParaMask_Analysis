#!/bin/bash
java -jar ~/scripts/ParaMask_simulation_scripts/ParaMask_Cluster_Seeds.jar \
	--cov GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03ES04Paramask.10PerMissing.snps.b.vcf.finalParaMaskIn.het.cov.stat.chr1.txt \
	--het ES0304_run_finalEMresults.chr1.het \
	--cutoff 254 \
	--covgw GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03ES04Paramask.10PerMissing.snps.b.vcf.finalParaMaskIn.het.cov.gw.txt \
	--range 1,29724344
