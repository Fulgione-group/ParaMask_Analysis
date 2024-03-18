install.packages("rtracklayer")
library(rtracklayer)

## import the bed file
bed.ranges <- import.bed("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Repeat_files/ES03_014_repeats_unpurged_assembly.bed")

## export as a gff3 file
export.gff3(bed.ranges,'ES03_014_repeat_regions.gff3')


ES04_bed.ranges <- import.bed("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Repeat_files/ES04_014_repeats_purged_assembly.bed")

## export as a gff3 file
export.gff3(bed.ranges,'ES04_014_repeat_regions.gff3')


### Liftover repeats with liftott ###
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load mambaforge/self-managed/v23.1.0-4
conda activate liftoff

GFF_file="/path/to/gff3" #ES03_014_repeat_regions.gff3
output_file="" #Repeats_ES03_on_reference.gff
intermediate_dir="" #intermediate_liftoff_files/
target_fasta="/path/to/ref.fa" #Reference_v5.1_only_chr.fasta
query_fasta="/path/to/assembly.fasta" #../../assemblies/ES03_014/final_assembly_minimap/ES03_014_assembly_only_chr.fasta

liftoff -g $GFF_file -o $output_file -dir $intermediate_dir $target_fasta $query_fasta
