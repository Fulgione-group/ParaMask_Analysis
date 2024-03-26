###
# Java programs 
###

Java programs to process the vcf, produce input files for analyses and compute summary statistics on the genomic data.

These java programs assume that you have the /java/ folder from GitHub in your home directory, and that the directory structure is the same as in /java/. 
If not so, you will need to change some paths (e.g. to the libraries in ./java/lib/)

To run these java programs, enter the ./java/projects directory and you will find:

### matrix_alpina.command
Reads a VCF file (download from EVA) and creates a SNP matrix filtered for coverage and quality. 

### sfs_lounch_2pops_alpina.command
Computes Joint and marginal allele frequency spectra for two populations, polarised to an outgroup (Arabis montbretiana in our case). 
Also outputs jsfs in dadi format.

in ./java/data you will find supporting files:


