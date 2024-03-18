#!/bin/bash
# sniffles
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load mambaforge/self-managed/v23.1.0-4
conda activate sniffles_env #sniffles v2.2
bsub -q multicore20 -R "rusage[mem=12000]" -M 15000 "/home/marimond/.conda/envs/sniffles/bin/sniffles --input ${ID}_reads_against_reference.bam --snf ${ID}.snf --threads 20 --vcf ${ID}_sniffles.vcf --reference ${reference}"

input_bam="/path/to/input/bam" #/biodata/dep_coupland/grp_fulgione/male/sv_calling/${ID}/${ID}_reads_against_referencev5.1.bam
ID="accession_ID"
sniffles --input input_bam --snf ${ID}_sniffles2_2_repeats.snf --threads 20 --vcf ${ID}_sniffles2.2_repeats.vcf --reference /biodata/dep_coupland/common/Arabis_alpina_resource/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta --threads 20 --minsvlen 35 --tandem-repeats /biodata/dep_coupland/grp_fulgione/male/sv_calling/Arabis_v5.1_repeats.bed --allow-overwrite
