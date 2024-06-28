#!/bin/bash
java -jar ParaMask_Cluster_Seeds.jar\
        --cov GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.goodformat.Q30LD5UD100K.vcf.het.cov.stat.chr9.txt\
        --het runs/PinkSalmon_nosex_EM2.7.1_EMresults_chr9.het\
        --covgw  GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.goodformat.Q30LD5UD100K.vcf.cov.gw.txt\
        --cutoff  569\
        --range 1,100203093
