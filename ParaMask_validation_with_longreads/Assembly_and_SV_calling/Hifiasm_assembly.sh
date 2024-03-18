#Hifiasm v0.19
#source /opt/share/software/scs/appStore/modules/init/profile.sh
#module load hifiasm/v0.19

table_file="/path/to/table_file/with/ID/and/raw_data_path.txt"
while IFS=$'\t' read -r col1 col2
do

ID="${col1}"
raw_data=${col2}"

#samples with heterozygosity > 0.1:

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/hifiasm_0.19
bsub -q multicore40 -R rusage[mem=60000] -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 ${raw_data}"

#samples with heterozygosity < 0.1:

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/hifiasm_0.19
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 -l0 ${raw_data}"

#after assembly transform gfa to fasta
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.fa



done < "$table_file"
