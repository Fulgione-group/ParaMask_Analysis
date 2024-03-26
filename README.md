# Overview
This repository encompasses all code that was used to test and analyze [ParaMask](https://github.com/Fulgione-group/ParaMask.git) results. For questions please contact us via [email](btjeng@mpipz.mpg.de).
<br>

## Simulations

Scripts to simulate a mosaic of single-copy and duplicated regions using forward genetic simulator with segmental duplications SeDuS ([Hartasánchez et al., 2015](https://academic.oup.com/bioinformatics/article/32/1/148/1742451)). These scripts were used to conduct simulations on a range of different duplicated sequence proportions, and with random mating and Inbreeding (F<sub>IS</sub>=0.9).

## SNP_calling
GATK ([DePristo et al., 2011](https://www.nature.com/articles/ng.806)) pipeline used to call SNPs for Spanish *A. alpina* populations that we tested ParaMask on.

## ParaMask_runs
Scripts and parameters used to run ParaMask  on Simulations and *A. alpina* populations.

## ParaMask_validation_with_long_reads
Scripts used to validate ParaMask calls on *A. alpina* populations based on long read sequences of two samples included in the Population. This mainly includes Assembly pipelines, Structural variant (SV) calling, Repeat calling, and scripts to analyze the overlap with ParaMask.    

## Jsfs_and_demographic_inference
Scripts to produce joint Site-Freqeuncy-Spectra (SFS), SFS based demographic inference using [dadi](https://dadi.readthedocs.io/en/latest/) ([Gutenkunst et al. 2009](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695)).

## Diversity_estimates
Scripts to estimate population mutation rate parameter $\theta$ and TajD on simulations and *A. alpina* populations.

## Fst_Fis_estimates
Scripts used to estimate F<sub>ST</sub> for *A. alpina populations* and F<sub>IS</sub>  for simulations and *A. alpina* populations.
