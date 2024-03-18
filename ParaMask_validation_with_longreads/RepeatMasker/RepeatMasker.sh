#!/bin/bash

### REPEATMASKER

fasta_assembly="/path/to/assembly/"
database_name="/database_name/"
#start building database
BuildDataBase -engine ncbi -name ${database_name} ${fasta_assembly}

#Model repeats and annotate repeats on assembly
RepeatModeler -pa 16 -engine ncbi -database ${database_name} 2>&1 | tee 00_repeatmodeler.log

#mask assembly with constructed database
RepeatMasker -gff -lib ${database_name}-families.fa ${fasta_assembly} > ${ID}_repeats_assembly.bed

awk '{print $5 "\t" $6 "\t" $7 "\t" $11}'  ES03_014_assembly_only_chr.fasta.out | sort | uniq | grep -v -E 'Unknown|Simple_repeat|Low_complexity|begin|position' > ES03_014_repeats.bed
awk '{print $5 "\t" $6 "\t" $7 "\t" $11}'  ES04_014_assembly_only_chr.fasta.out | sort | uniq | grep -v -E 'Unknown|Simple_repeat|Low_complexity|begin|position' > ES04_014_repeats.bed

