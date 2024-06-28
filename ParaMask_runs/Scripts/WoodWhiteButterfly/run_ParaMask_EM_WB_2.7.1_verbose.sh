#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate Paramask_env

#ulimit -v 50000000
#export OPENBLAS_NUM_THREADS=80 OMP_NUM_THREADS=80 MKL_NUM_THREADS=40
Rscript --vanilla ParaMask_EM_v0.2.7.1.R\
      	--het WB_allchrom.Q30LD5UD100K.vcf.het.stat.txt\
        --missingness 0.1\
	--outdir $DIR\
	--verbose\
	--ID WB_EM_2.7.1
