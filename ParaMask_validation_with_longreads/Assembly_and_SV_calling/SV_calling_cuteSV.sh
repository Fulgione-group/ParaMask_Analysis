# cuteSV
bsub -q multicore20 -R "rusage[mem=22000]" -M 25000 "/home/marimond/.conda/envs/bases/bin/cuteSV ${col1}_reads_against_reference.bam ${reference} ${ID}_cutesv.vcf ${work_dir} -l 35"
