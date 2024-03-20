#!/bin/bash

#Liftover repeats of ES03 and ES04 to reference to fit coordinates (using liftoff, input created gff file of repeats)

#source /opt/share/software/scs/appStore/modules/init/profile.sh
#module load mambaforge/self-managed/v23.1.0-4
conda activate liftoff

GFF_file="/path/to/gff3" #ES03_014_repeat_regions.gff3
output_file="/Path/to/output_dir/" #Repeats_ES03_on_reference.gff
intermediate_dir="/Path/to/intermediate_files/" #intermediate_liftoff_files/
target_fasta="/path/to/ref.fa" 
query_fasta="/path/to/assembly.fasta" 
liftoff -g $GFF_file -o $output_file -dir $intermediate_dir $target_fasta $query_fasta
