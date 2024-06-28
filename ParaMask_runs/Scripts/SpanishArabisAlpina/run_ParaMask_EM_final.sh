#!/bin/bash
source /netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/anaconda3/etc/profile.d/conda.sh
conda activate Paramask_env

ulimit -v 50000000
export OPENBLAS_NUM_THREADS=40 OMP_NUM_THREADS=40 MKL_NUM_THREADS=40
Rscript --vanilla ./ParaMask_EM_v2.4.R "--hetfile" "GATK4.2_1000Genomes_chrall.filteredQ30LD5UD100K.final.ES03ES04Paramask.10PerMissing.snps.b.vcf.finalParaMaskIn.het.het.stat.txt" \
"--outpath" "$DIR" "--missingness" "0.1" "--verbose" "--ID" "ES0304_run_final" "--dist_em_rep" "100"
