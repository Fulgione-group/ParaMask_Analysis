# Overview

We simulated duplications interdispersed among single-copy regions for different proportions of duplications $p$ from 0% to 50%. For this we first sampled sequence lengths of duplication from a exponential distribution and single-copy sequence lengths from a geometric distribution for a total sequence length of 1 Mbp. Than we ran simulations for each duplicated and single-copy region seperately using segmental duplication simulator [SeDuS](https://academic.oup.com/bioinformatics/article/32/1/148/1742451). After simulations we collapsed duplicated regions by converting sites with a derived alleles on exactly one of the copies to heterozygotes. Next we simulated coverage for each site and genotype by sampling from a negative binomial with overdispersion parameter $\theta$<sub>NB</sub> = 8, and mean 10 or 20 for single-copy regions and duplications respectively. Sites with simulated coverage less than 3 were excluded from the analysis. Lastly we simulate read ratios by sampling from a binomial with proportion 0.5 and a binomial population size equal to the simulated coverage. 


## Sample sequence length
### Script: Sample_duplication_length_cmd.R
This R script will sample sequence lengths of duplication from a exponential distribution with a mean of 1000 bp and single-copy sequence lengths from a geometric distribution with mean = (L<sub>total</sub>-p)/(N<sub>duplications</sub> + 1) for a total sequence length L<sub>total</sub> of 1 Mbp. 
