#!/bin/bash

awk -v subsample="$1" 'BEGIN {srand(); for (i=1; i<=100; i++) a[i]=i; for (i=100; i>1; i--) { j=int(rand()*i)+1; t=a[i]; a[i]=a[j]; a[j]=t; } for (i=1; i<=subsample; i++) {b[a[i]+9]=1};}{if($0!~/^##/){printf"%s", $1;for(i=2;i<=9;i++){printf"\t%s", $i};for(i=10;i<=109;i++){if(b[i]==1){printf"\t%s", $i}}; printf"\n"}else{print $0}}' $2 > $(echo $2 |  cut -d'.' -f1)"_sub"$1".vcf"
