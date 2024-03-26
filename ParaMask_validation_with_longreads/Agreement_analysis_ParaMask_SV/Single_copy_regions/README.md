The script _Filter_for_seeds_file_preparation.sh_ is employed to filter all ParaMask cluster name files, retaining multi-copy cluster names identified in ES03-014 and/or ES04-014 exclusively.

Subsequently the script _Filter_for_ES03_ES04_seeds_only.R_ selectively filters ParaMask calls using the filtered cluster names, identifying only multicopy SNPs of ES03-014 and/or ES04-014, as well as multicopy clusters with a minimum length of 35 base pairs.

To facilitate the calculation of agreement between ParaMask and long-read SV callers within single-copy regions, the script _Function_agreement_single_copy.R_ is utilized.

This function is applied to all regions that are not annotated as structural variants by either cuteSV or sniffles. These regions are generated using the _script Generate_single_copy_regions.R_.
