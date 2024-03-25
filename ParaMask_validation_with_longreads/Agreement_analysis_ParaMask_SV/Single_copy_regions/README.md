The script Filter_for_ES03_ES04_seeds_only.R filters the ParaMask calls for multicopy SNPs called for ES03-014 and/or ES04-014 as well as for cluster with a minimum of 35bp only.
With the script Filter_for_seeds_file_preparation.sh all ParaMask cluster files are filtered for ES03-014 and ES04-014 seeds only.

To define a function to calculate the agreement between ParaMask and long-read SV callers in single-copy regions the script Function_agreement_single_copy.R is used.

With the script Generate_single_copy_regions.R, single-copy regions based on all merged not overlapping structural variants are generated. 
