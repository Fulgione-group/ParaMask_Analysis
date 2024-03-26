#!/bin/bash
acr=$(awk -v l=$1 'BEGIN{print (10*l/5000)}')

/netscratch/dep_coupland/grp_fulgione/bastiaan/software/SeDuS_1.10_exe Sim_SV_b$2_l$1_c0.05_rep$3\
 -s 1\
 -z 100\
 -k 1000\
 -n 5250\
 -b $1\
 -u 0.000525\
 -r $acr\
 -c 0.05
