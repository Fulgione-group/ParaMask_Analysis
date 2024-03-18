#!/bin/bash

BAM_file ="/Path/to/Bam_file.bam"
ID="Sample_ID"

samtools depth -a ${BAM_file} > ${ID}_coverage.coverage
