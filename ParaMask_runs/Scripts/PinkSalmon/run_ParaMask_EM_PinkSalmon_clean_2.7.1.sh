#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate Paramask_env

#ulimit -v 50000000
#export OPENBLAS_NUM_THREADS=80 OMP_NUM_THREADS=80 MKL_NUM_THREADS=40
Rscript --vanilla /home/btjeng/scripts/ParaMask_beta_release/ParaMask/ParaMask/ParaMask_EM_v0.2.7.1.R\
      	--het GL_PinkSalmon_even_GATKfilter.GATKrem.snps.biallelic.goodformat.Q30LD5UD100K.vcf.het.stat.txt\
        --missingness 0.1\
	--outdir Validation_ParaMask/Pink_salmon\
	--ID PinkSalmon_nosex_EM2.7.1
