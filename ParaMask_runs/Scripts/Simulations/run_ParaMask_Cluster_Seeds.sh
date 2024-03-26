#!/bin/bash
java -jar ParaMask_Cluster_Seeds.jar \
	--cov $1 \
	--het $2 \
	--cutoff $3 \
	--covgw $4 \
	--range 1,1000000
