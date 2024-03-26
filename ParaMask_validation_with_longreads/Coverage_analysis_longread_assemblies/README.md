This directory contains two scripts for a coverage analysis of whole genome assemblies generated with long-read data.

Utilizing a mapping file in BAM format, the script _Calculate_coverage_per_base.sh_ determines the depth per base and subsequently computes the coverage within a 1500-base pair window.

To annotate potential coverage peaks and drops, _Annotation_of_coverage_peaks.R_ is used. This script additionally needs results of the RepeatMasker pipeline as input. 
