# Overview

We simulated duplications interdispersed among single-copy regions for different proportions of duplications $p$ from 0% to 50%. For this we first sampled sequence lengths of duplication from an exponential distribution and single-copy sequence lengths from a geometric distribution for a total sequence length of 1 Mbp. Than we ran simulations for each duplicated and single-copy region seperately using segmental duplication simulator [SeDuS](https://academic.oup.com/bioinformatics/article/32/1/148/1742451). After simulations we collapsed duplicated regions by converting sites with a derived alleles on exactly one of the copies to heterozygotes. Next we simulated coverage for each site and genotype by sampling from a negative binomial with overdispersion parameter $\theta$<sub>NB</sub> = 8, and mean 10 or 20 for single-copy regions and duplications respectively. Sites with simulated coverage less than 3 were excluded from the analysis. Lastly we simulate read ratios by sampling from a binomial with proportion 0.5 and a binomial population size equal to the simulated coverage. 
<br />
<br />
<br />
## Sample sequence length
### Script: Sample_duplication_length_cmd.R
This R script samples sequence lengths of duplications from an exponential distribution with a mean of 1000 bp for a total duplicated sequence length of p * 1 Mbp, where p is the proportion of duplications, and samples single-copy sequence lengths from a geometric distribution with mean = (1 Mbp * (1-p))/(N<sub>duplications</sub> + 1), where N<sub>duplications</sub> is the number of duplications drawn until the total duplicated sequence length is reached. The script outputs a text table with single-copy and multicopy regions respectively.
<br />
<br />
<br />
## Simulation using SeDuS

The full steps for simulations with random mating and inbreeding are summarized in the the Pipeline_sim.sh and Pipeline_sim_withSelfing.sh scripts respectively. The most important steps are explained here:
<br />
### SeDus_simSC_var.sh and SeDus_simSV_var.sh
These scripts dynamically run SeDuS with a sequence length as input. SeDus then performs forward simulations.

```bash
#!/bin/bash
acr=$(awk -v l=$1 'BEGIN{print (100*l/5000)}')

/netscratch/dep_coupland/grp_fulgione/bastiaan/software/SeDuS_1.10_exe Sim_SC_b$2_l$1_c0.05_rep${3}\
 -s 1\
 -z 100\
 -k 1000\
 -n 10000\
 -b $1\
 -u 0.001\
 -r ${acr}\
 -c 0.05
```
Where "-b" is the sequence length and "-r" is the recombination rate adjusted for the sequence length, "-u" is the fixed population mutation rate parameter $\theta$, "-n" the simulated population size, and "-c" is the population-scaled interlocus gene conversion rate.

For simulations with Inbreeding (F<sub>IS</sub>=0.9) we adjusted population mutation rate parameter theta and recombination rate according to [Nordborg 2000](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1460950/)

For single-copy sequences only block 1 outputed by SeDuS is processed and for duplication blocks 0 and 2.
<br />
### genotype_to_mat_SC_new.sh  genotype_to_mat_SV2_new.sh
These awk scripts store the original genotype information. Homozygotes, heterozygotes, homozygotes at copies of duplication, and heterozygotes at copies of duplications are encoded differently.  
<br />
### concat_SV_mat.sh and concat_SC_SV.sh
These scripts collapse both copies of the duplications and concatinate the single-copy sequences and duplications.
<br />
<br />
<br />
## Simulation of coverage and read ratiosS
### Simulate_coverage_allele_depth
This R script simulates coverage for each site and genotype by sampling from a negative binomial with overdispersion parameter $\theta$<sub>NB</sub> = 8, and mean 10 or 20 for single-copy regions and duplications respectively. Sites with simulated coverage less than 3 were excluded from the analysis. We simulate the number of reference alleles for heterozygote genotypes by sampling from a binomial distribution with a binomial population size equal to the simulated coverage. For heterozygote genotypes at single-copy regions and for homozygote differences between duplicated copies, the binomial probability to observe the derived allele was 0.5. For sites that are heterozygote on one duplicated copy and respectively homozygote ancestral or derived on the other, the probability was 0.25 or 0.75.


