#!/bin/bash
java -jar PrepareParaMaskInput_fromVCF.jar \
-VF \
--vcf $1 \
--missingness 0.1

