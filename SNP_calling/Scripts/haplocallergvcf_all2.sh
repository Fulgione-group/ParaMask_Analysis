#!/bin/bash/

dir=$1
echo "dir"
cd $dir
samples=$2
bootrun=$3
maxjob=$(expr $(expr $4 - 1 ) / 8 )

gatk2='/netscratch/dep_coupland/grp_coupland/bioinformatics/bastiaan/software/gatk-4.2.0.0/gatk'

haplocaller () {
	local j="$1"
	samtools index $j &&
	Out=$(echo $j | rev | cut -d '_' -f2- | rev)
	$gatk2 --java-options "-XX:ConcGCThreads=2 -XX:ActiveProcessorCount=8 -Xmx10g" HaplotypeCaller \
   		-R /netscratch/dep_coupland/grp_fulgione/bastiaan/data/reference/Alpina_V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta \
    		-I $j \
    		-O ${Out}_b${bootrun}_allsites.g.vcf \
    		-ERC GVCF \
		--output-mode EMIT_ALL_ACTIVE_SITES \
    		--native-pair-hmm-threads 8
}


for i in $(eval ls ${samples}_mergebamalignment_dm.bam)
do
        haplocaller $i &
        background=( $(jobs -p) )
        if (( ${#background[@]} == $maxjob )); then
                        wait -n
        fi
done
wait

