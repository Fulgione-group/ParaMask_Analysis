#!/bin/bash
cd $3
awk 'BEGIN{FS=OFS="\t"}{if(NR==FNR){if(FNR>1){p[$4]=($1-1)}}else{$1=($1+p[$NF]); print $0}}' $1 $2> $(echo $2 | rev | cut -d'.' -f2- |rev)_cpos.txt
