#!/bin/bash

##Error rate for heterozygotes is 0.01 either homozygous
awk 'BEGIN{
	FS=" ";
	OFS="\t";
}{
	if($1!=0){
		split(FILENAME, z,"_");
		block=z[6];
		for(i=2;i<=NF;i++){
			if(i%2==0){
				a=$i
			}else{
				b=$i;
				#het
				if((a+b)==1){
					g[((i-1)/2)]=1
				}else if((a+b)==2){
					g[((i-1)/2)]=2
				}else{
					g[((i-1)/2)]=0
				}
				b=0;
				a=0
			}
		};
		printf"%s\t", $1
		for(i=1; i<=(length(g)-1); i++){
			printf"%s\t", g[i];
		}
		printf"%s\t%s\t%s\n",g[length(g)], (NF-1), "SC_"block;
	}
}' $1 > $2

