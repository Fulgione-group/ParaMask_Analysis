#!/bin/bash

BAM_file ="/Path/to/Bam_file.bam"
ID="Sample_ID"
WINDOW_SIZE=1500
OUTPUT_BED="$SCRIPT_DIR/mean_coverage.bed"

samtools depth -a ${ID}_pb.paf.gz | \
    awk -v window="$WINDOW_SIZE" 'BEGIN{OFS="\t"; start=1; end=window} {
        while (start < $2) {
            total=0; count=0
            for (i=start; i<=end; i++) {
                total+=$3; count++
            }
            if (count > 0) {
                print $1, start, end, total/count
            }
            start+=window; end+=window
        }
    }' > "$OUTPUT_BED"

