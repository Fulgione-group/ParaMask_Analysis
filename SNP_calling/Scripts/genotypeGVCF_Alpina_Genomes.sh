#!/bin/bash

dir=$1
echo "dir"
cd $dir
samples=$2
gatk2='/netscratch/dep_coupland/grp_fulgione/bastiaan/software/gatk-4.2.1.0/gatk'


$gatk2 --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -XX:ConcGCThreads=2 -XX:ActiveProcessorCount=4 -Xms10G -Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
    -R /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
    -V gendb:///netscratch/dep_coupland/grp_fulgione/bastiaan/data/GATK_DBs/AlpinaGenomes \
    --use-new-qual-calculator \
    -O GATK4.2_alpina_15_02_2023_chr1.vcf \
    -all-sites\
    --tmp-dir /scratch\
    -L chr1
