# Overview 
We used [GATK](https://www.nature.com/articles/ng.806) [v4.2.0.0](https://gatk.broadinstitute.org/hc/en-us/sections/360012354372-4-2-0-0) for SNP calling in the two *Arabis alpina* populations. Read groups were assigned, adapters were soft clipped, duplicates were marked and reads were aligned to the pajares reference version [v5.1](http://www.arabis-alpina.org/data/ArabisAlpina/assemblies/V5.1/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta.gz) ([W.-B. Jiao and K. Schneeberger, 2017](https://www.sciencedirect.com/science/article/pii/S1369526616301315?via%3Dihub)) ) using the mem algorithm of bwa v0.7.17 ([H. Li and R. Durbin, 2009](https://pubmed.ncbi.nlm.nih.gov/19451168/))

<br>
<br>
<br>
# Pipeline:
## fastqtosam_all.sh

Script convert fastq files to unmapped sam files and assign read groups from header information (valid for DNBseq headers 2020). 
<br>
## markadapters_all.sh and discount_all.sh
\*Note: Adapters were already cleaned by the sequencing provider, Beijing Genomics Institute (BGI). 
Scripts to markadapters and soft clip them, supplied with five prime and three prime adapters from DNBseq.
<br>
## bwamem_all.sh
Script to map samples to the *A. alpina* pajares reference.
<br>
## mergebamalignment_all.sh and sortmabam_all.sh
Script to merge mapped bam file and raw unmapped sam files and sort the merged bam file
<br>
## markadapters
