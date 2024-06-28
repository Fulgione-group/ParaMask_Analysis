#!/bin/bash
java -jar PrepareParaMaskInput_fromVCF.jar\
        --vcf WB_allchrom.Q30LD5UD100K.vcf.gz\
        --missingness 0.1

