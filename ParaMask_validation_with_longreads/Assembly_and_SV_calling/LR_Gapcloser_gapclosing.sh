#!/bin/bash
export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag_nucmer/ragtag_wo_corr/ragtag.scaffold.fasta -l SE02_002_raw.fasta -o LR_Gapcloser/"

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
lsfenv
while read -r col1 col2 _; do
  samtools faidx /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta
  echo "Index of ${col1} created"
  seqtk seq -a ${col2} > /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/${col1}_raw.fasta
  export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
  bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta -l /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/${col1}_raw.fasta -m 1000 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser"
done < "$input_file"
