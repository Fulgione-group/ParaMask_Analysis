# Overview

For each site the number of genotypes are randomly projected down to 90% in order to account for missingsites. Then we calculated diversity estimators $\theta$<sub>W</sub> and $\theta$<sub>$\pi$</sub>, and Tajimas' D. $\theta$<sub>W</sub> and $\theta$<sub>$\pi$</sub> are scaled by the density of sites.

# Scripts

Scripts are subdivided into two folders "Simulations" and "SpanishArabisAlpina" for simulated and genomic data for *A. alpina* populations ES03 and ES04 respectively. For the two *A. alpina* populations diversity estimates $\theta$<sub>W</sub> and $\theta$<sub>$\pi$</sub> were scaled to number of sites retained after filtering all sites (including invariant) for a maximum missingness of 10% (can be found in the scripts).

# Bytecode
Compiled java bytecode for downsampling and calculating diversity estimates.

# src
Java source code  for downsampling and calculating diversity estimates.
