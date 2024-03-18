#!/bin/bash
reference="/netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta"
assembly="ASSEMBLY"
ID=""
mapper="PATH_TO_MINIMAP_OR_NUCMER"

export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
# without chimeric contig correction
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold ${reference} ${assembly} -o ${ID}_RagTag -C --aligner ${mapper} -t 20"

export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/merged_flye_hifiasm/merged_out.fasta -C -t 20"
