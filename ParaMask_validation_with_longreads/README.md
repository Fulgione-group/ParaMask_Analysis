# ParaMask validation with longreads
This directory contains scripts for validating the output of ParaMask against PacBio Hifi long-read structural variant (SV) calls for two accessions belonging to the respective Arabis alpina population.

## Assembly and SV calling
The directory assembly and SV calling contains all scripts used to assemble PacBio Hifi longreads of ES03 and ES04 into whole genome sequences and to perform SV calling. 

## Coverage Analysis longread assemblies 
The directory Coverage Analysis longread assemblies contains all scripts used to annotate high coverage peaks inferred based on the alignment between assembled scaffolds with their corresponding longreads. Additionally, files of the RepeatMasker Analysis and the results of the Blast Search of the extracted high coverage peaks were used as input for the annotation. 

## Fit coordinate systems

## Agreement analysis ParaMask SV
The directory Agreement analysis contains all scripts to validate the overlap between ParaMask and SV calls seperated for single-copy and multicopy regions. Each subdirectory contains a R-script to define a function how to calculate the agreement in the particular region. The single-copy directory additionally contains R-scripts to filter the ParaMask calls for ES03 and ES04 seed calls only. 

## RepeatMasker

## Agreement analysis ParaMask SV 

## Repeat SV overlap analysis

