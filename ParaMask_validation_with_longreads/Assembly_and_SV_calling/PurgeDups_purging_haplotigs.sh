#!/bin/bash
# Purging of haplotigs with purge_dups
# Heterozygosity threshold: 0.1
raw="path/to/data"
ID="ID"
assembly="assembly.fasta"
tool="flye"

#Map reads to assembly and calculate depth +
minimap2 -xasm20 ${assembly} ${raw}  | gzip -c - > ${col1}_${tool}_pb.paf.gz
pbcstat ${ID}_${tool}_pb.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

#split assembly and do self-alignment
split_fa ${assembly} > ${assembly}.split
minimap2 -xasm5 -DP  ${assembly}.split  ${assembly}.split | gzip -c - >  ${assembly}.split.self.paf.gz

#purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${assembly}.split.self.paf.gz > dups.bed 2> purge_dups.log

#get purged sequences and haplotigs
get_seqs -dups.bed ${assembly}


#bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bash purge_dups.sh"
