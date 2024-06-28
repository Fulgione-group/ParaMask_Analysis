#!/bin/bash
java -jar ParaMask_Cluster_Seeds.jar\
	--cov WB_allchrom.Q30LD5UD100K.chr20.vcf.het.cov.stat.txt\
	--het WB_EM_2.7.1_EMresults.chr20.het\
	--covgw WB_allchrom.Q30LD5UD100K.vcf.cov.gw.txt\
	--cutoff  300\
	--range 1,16307044
