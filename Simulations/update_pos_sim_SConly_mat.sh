#!/bin/bash
awk 'BEGIN{FS=OFS="\t";c=0}{if(FNR==NR){l[$2]=c; c+=$1}else{$1+=l[$NF];print $0}}' $1
