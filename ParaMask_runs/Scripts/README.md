# Overview 
Scripts to run Paramask on Simulations and Spanish *Arabis alpina* populations. For all runs we used biallic SNPs with genotype quality >=30 and min depth >=5.
<br>
## Simulations
for simulations we used default settings with maximum 10% missingness per site for all ParaMask steps.
## SpanishArabisAlpina, WoodWhiteButterfly and  PinkSalmon
For the *A. alpina* population genomic data we used default setting for all parameters with maximum missingness of 10% per site in all Paramask steps, except that we set the number of iterations in the distance EM algorithm to 100.  
For WoodWhiteButterfly population genomic data we used default setting for all parameters with maximum missingness of 10% per site in all Paramask steps
For PinkSalmon population genomic data we used default setting for all parameters with maximum missingness of 10% per site in all Paramask steps. The script for painting chromosomes according to the copynumber status is included.
For Chinook Salmon we used a maximum missingness of 30% comparable to inference by HDplot and rCNV. Clustering step for this RADseq data was ommitted.
## Evaluation
Scripts to evaluate the Overlap between simulated multicopy regions and multicopy regions identified by ParaMask
Scripts to analyze the ParaMask results for genomic data of Arabis alpina, Pink Salmon, and Wood White Butterfly
Script to evaluate the recall of known multicopy SNPs in Chinook Salmon
