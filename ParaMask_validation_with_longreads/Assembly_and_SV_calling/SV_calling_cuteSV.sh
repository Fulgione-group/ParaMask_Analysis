#/bin/bash

# SV calling using cuteSV

CuteSV="/Path/to/cuteSV/"
ID="Accession_ID"
reference="/Path/to/reference.fasta"
work_dir="/Path/to/working_directory/"

${CuteSV}/bin/cuteSV ${ID}_reads_against_reference.bam ${reference} ${ID}_cutesv.vcf ${work_dir} -l 35"
