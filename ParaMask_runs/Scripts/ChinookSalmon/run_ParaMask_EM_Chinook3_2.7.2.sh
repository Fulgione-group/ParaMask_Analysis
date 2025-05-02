#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate Paramask_env

#ulimit -v 50000000
#export OPENBLAS_NUM_THREADS=80 OMP_NUM_THREADS=80 MKL_NUM_THREADS=40
Rscript --vanilla ParaMask_EM_v0.2.7.2.R\
        --het Chinook_sequenceReads_IDCHR_v3.vcf.het.stat.txt\
        --missingness 0.3\
        --outdir Validation_ParaMask/Chinook/\
        --ID ParaMask_Chinook3

