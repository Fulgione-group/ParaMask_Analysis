#!/bin/bash
java -jar ~/scripts/ParaMask_simulation_scripts/PrepareParaMaskInput_fromVCF.jar \
-VF \
--vcf $1 \
--missingness 0.1

