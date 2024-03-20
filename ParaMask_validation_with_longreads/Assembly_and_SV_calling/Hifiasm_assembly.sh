#!/bin/bash

#Hifiasm v0.19
#source /opt/share/software/scs/appStore/modules/init/profile.sh
#module load hifiasm/v0.19

table_file="/path/to/table_file.txt/"   #including accession IDs (col 1) and paths to raw reads (col 2)
  while IFS=$'\t' read -r col1 col2
    do

    ID="${col1}"
    raw_data=${col2}"
    working_dir="/path/to/working_dir/"

#samples with heterozygosity > 0.1: (purging)

cd ${working_dir}
bsub -q multicore40 -R rusage[mem=60000] -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 ${raw_data}"

#samples with heterozygosity < 0.1: (no purging -l0)

cd ${working_dir}
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 -l0 ${raw_data}"

#after assembly transform gfa to fasta
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.fa



done < "$table_file"
