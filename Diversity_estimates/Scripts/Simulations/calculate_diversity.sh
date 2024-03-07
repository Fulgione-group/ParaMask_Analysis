#!/bin/bash
java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep1_mat_cpos.vcf\
        --cutoff 0.1\
        --seqlength 1000000\
        --id 0.1_Collapsed\
        --replicate 1

java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep2_mat_cpos.vcf\
        --cutoff 0.1\
        --seqlength 1000000\
        --id 0.1_Collapsed\
        --replicate 2


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep3_mat_cpos.vcf\
        --cutoff 0.1\
        --seqlength 1000000\
        --id 0.1_Collapsed\
        --replicate 3

java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep1_mat_cpos_sc.vcf\
        --cutoff 0.1\
        --seqlength 898353\
        --id 0.1_SC\
        --replicate 1


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep2_mat_cpos_sc.vcf\
        --cutoff 0.1\
        --seqlength 896599\
        --id 0.1_SC\
        --replicate 2


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep3_mat_cpos_sc.vcf\
        --cutoff 0.1\
        --seqlength 899665\
        --id 0.1_SC\
        --replicate 3


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep1_mat_cpos_ParaOut.vcf\
        --cutoff 0.1\
        --seqlength 889912\
        --id 0.1_ParaOut\
        --replicate 1


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep2_mat_cpos_ParaOut.vcf\
        --cutoff 0.1\
        --seqlength 887761\
        --id 0.1_ParaOut\
        --replicate 2


java -jar /home/btjeng/Data/SeDuS/local_SeDuS_runs/ParaMask_Calculate_Diversity.jar \
        --vcf Simulations_SV_SC_rep3_mat_cpos_ParaOut.vcf\
        --cutoff 0.1\
        --seqlength 890929\
        --id 0.1_ParaOut\
        --replicate 3


