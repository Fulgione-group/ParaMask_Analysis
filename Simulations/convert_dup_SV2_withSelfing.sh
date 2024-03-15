#!/bin/bash

#find duplicated position because of gene conversion
awk -v Fis=$4 'BEGIN{empty=1;
	probhet=(1-Fis)
        probhom=(1-probhet)/2
        FS=" ";
        OFS="\t";
        srand(systime() + int(rand() * 1000000) + PROCINFO["pid"]);
}{
	if(FNR==NR){
		if($1!=-1 && $1!=0){
			p[FNR]=$1;
			l[FNR]=$0
		}
	}else{
		split(FILENAME, z,"_");
        	block=z[6];
		for(i=1;i<=length(p); i++){
			if($1==p[i]){
				srand();
				empty=0;
				split(l[i],la," ");
				#print "yeah";
				printf"%s\t", $1;
				for(j=2;j<=NF;j++){
					#sum over all 4 genotypes
					gsum1=$j+$(j+1)
					gsum2=la[j]+ la[(j+1)];
					#transform because of selfing
					if(gsum1==1){
						r=rand();
                                       		if(r<probhet){
                                        		gsum1=1
                                       		}else if((r-probehet)<probhom){
                                        		gsum1=2
                                        	}else {
                                        		gsum1=0
                                       		}
					}
                                        if(gsum2==1){
                                                r=rand();
                                                if(r<probhet){
                                                        gsum2=1
                                                }else if((r-probehet)<probhom){
                                                        gsum2=2
                                                }else {
                                                        gsum2=0
                                                }
                                        }
					gsum=gsum1+gsum2
					#ref call
					if(gsum==0){
						g=0
					#alt call
					}else if(gsum==4){
						g=2
					#clear het , 2% chance beeing fixed, 98%chance beeing het
					}else if(gsum==2){
						g=1
					#sometimes ref (40%) or het (40%) or Na (20 %)
					}else if(gsum==1){
						g=3
					#sometimes alt (40%) or het (40%) or Na (20%)
					}else if(gsum==3){
						g=4
					}
					;printf"%s\t",g;
					j++
				}
			printf"%s\t%s\n",(NF-1), "SV_"block
			}
		}
	}
}END{
	if(empty==1){
		printf"%s\t%s\n",-1, (NF-1),"SV_"block;
	}
}' $1 $2 > $3
