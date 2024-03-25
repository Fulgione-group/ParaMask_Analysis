#!/bin/bash

##Error rate for heterozygotes is 0.01 either homozygous
awk 'BEGIN{
	probhet=SR/(2 - SR)
        probhom=(1-probhet)/2
	FS=" ";
	OFS="\t";
	srand(systime() + int(rand() * 1000000) + PROCINFO["pid"]);
}{
	if(FNR==NR){
		if($1!=-1){
			pos[FNR]=$1
		}
	}
	else{
		if($1!=0){
			duppos=0;
			for(i=1; i<= length(pos); i++){
				if(pos[i]==$1){
					duppos=1;
				}
			}
			if(duppos==0){
				split(FILENAME, z,"_");
				block=z[6]
				for(i=2;i<=NF;i++){
					if(i%2==0){
						a=$i
					}else{
						b=$i;
						#if heterozygous 0.4 prob to be collapsed heterozygous, 0.4 prob to be homozygous ancestral and 0.2 prob to be missing 
						if((a+b)==1){
                        	                	g[((i-1)/2)]=3
						#homozygous derived 0.98 prob to appear heterozygous, 0.01 prob to appear homozygous ancestral or derived 
						}else if((a+b)==2){
	                                	       	g[((i-1)/2)]=1
						}else{
                        	       		        g[((i-1)/2)]=0
						}
						b=0;
						a=0;
					}
				};
        			printf"%s\t", $1
        			for(i=1; i<=(length(g)-1); i++){
        			        printf"%s\t", g[i];
        			}
        			printf"%s\t%s\t%s\n",g[length(g)], (NF-1), "SV_"block;
			}
		}
	}
}' $1 $2 > $3
