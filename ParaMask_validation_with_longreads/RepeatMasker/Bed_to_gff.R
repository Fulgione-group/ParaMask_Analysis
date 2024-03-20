install.packages("rtracklayer")
library(rtracklayer)

## import the bed file
bed.ranges <- import.bed("~/Path/to/Repeat_files/ES03_014_repeats_unpurged_assembly.bed")

## export as a gff3 file
export.gff3(bed.ranges,'ES03_014_repeat_regions.gff3')


ES04_bed.ranges <- import.bed("~/Path/to/Repeat_files/ES04_014_repeats_purged_assembly.bed")

## export as a gff3 file
export.gff3(bed.ranges,'ES04_014_repeat_regions.gff3')



