#!/bin/bash
java -jar ~/scripts/ParaMask_simulation_scripts/ParaMask_Cluster_Seeds.jar \
	--cov $1 \
	--het $2 \
	--cutoff $3 \
	--covgw $4 \
	--range 1,1000000
