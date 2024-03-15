#!/bin/bash


###
#		VCF combined to SNP matrix
###
vcf=/srv/netscratch/dep_coupland/grp_fulgione/andrea/alpina_2023/GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.b.vcf.gz
# bgzip -cd ${vcf}.b.gz > ${vcf}
matrix=${vcf}_e3e4_paramask_superMatrix.txt
vcfNew=${vcf}.b.gz_e3e4_paramask.b.vcf

java -Xmx4G -classpath ~/java/lib/junit.jar:~/java/lib/jbzip2-0.9.1.jar:~/java/lib/args4j-2.0.12.jar:~/java/lib/commons-compress-1.0.jar:~/java/lib/gson-1.6-javadoc.jar:~/java/lib/gson-1.6-sources.jar:~/java/lib/gson-1.6.jar:bin c.e.data_processing.VcfCombined_to_snpMatrix_alpina ${vcf} ${1} ${matrix}

java -Xmx4G -classpath ~/java/lib/junit.jar:~/java/lib/jbzip2-0.9.1.jar:~/java/lib/args4j-2.0.12.jar:~/java/lib/commons-compress-1.0.jar:~/java/lib/gson-1.6-javadoc.jar:~/java/lib/gson-1.6-sources.jar:~/java/lib/gson-1.6.jar:bin c.e.data_processing.Matrix_snpsOnly ${matrix}_c5q25_chr${1}.txt


