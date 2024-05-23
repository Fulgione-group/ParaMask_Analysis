#!/bin/bash
bam_list="/PATH/TO/BAMLIST.txt"
snp_position_file="/PATH/TO/SNP_to_TEST.txt"
ngs_in_file="/PATH/TO/NGS_FILE.pileup"

# Create pileup file 
samtools mpileup -b ${bam_list} -l ${snp_position_file} -q 0 -Q 0--ff UNMAP,DUP  > ${ngs_in_file}

# Split pileup file per chromosome
chromosomes=$(awk '{print $1}' "$ngs_in_file" | sort | uniq)
for chr in $chromosomes; do
    output_file="${chr}.pileup"
    awk -v chr="$chr" '$1 == chr' "$input_file" > "$output_file"
    echo "Created file: $output_file"
done

# Rung ngsparalog calcLR function per chromosome 
for chr in $chromosomes; do
  export PATH=$PATH:/netscratch/dep_coupland/grp_fulgione/male/ngsParalog
  current_chr_file="${chr}.pileup"
  outfile="/netscratch/dep_coupland/grp_fulgione/male/ParaMask/ngsParalog/NGS_out_${chr}.lr"
  ngsParalog calcLR -infile ${current_chr_file} -outfile ${outfile} -minQ 20 -minind 25 -mincov 1
done
