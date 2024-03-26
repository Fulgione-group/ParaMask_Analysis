#!/bin/bash
source PATH_TO_CONDA/conda.sh
conda activate Paramask_env #for conda environment with all prerequisites


##$1= hetfile, $2=output dir, $3=ID

ulimit -v 50000000
export OPENBLAS_NUM_THREADS=40 OMP_NUM_THREADS=40 MKL_NUM_THREADS=40
Rscript --vanilla PATH_TO_PARAMASK/ParaMask_EM_v2.4.R "--hetfile" $1 \
"--outpath" $2 "--missingness" "0.1" "--verbose" "--ID" $3
