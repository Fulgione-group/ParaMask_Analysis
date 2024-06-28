#!/bin/bash
java -jar PrepareParaMaskInput_fromVCF.jar\
        --vcf GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.goodformat.Q30LD5UD100K.vcf.gz\
        --missingness 0.1

