#!/bin/bash
#Filter for ES03 / ES04 clusters only: 

# Function to process clusters for a given chromosome and seed
process_clusters() {
    awk -v chr="$1" -v seed="$2" '
        BEGIN { n=0; l=0; t=0; c=0 }
        {
            if (FNR > 1) {
                if ($3 != l) {
                    if (t == 1) {
                        c += n;
                        for (i=1; i<=length(k); i++) {
                            print k[i];
                        }
                    }
                    t=0; n=0; split("", k);
                }
                l=$3; n++; k[n]=$0;
                split($NF, a, ":");
                for (i=1; i<=length(a); i++) {
                    if (a[i] == seed) {
                        t=1;
                    }
                }
            }
        }
        END {
            if (t == 1) {
                c += n;
                for (i=1; i<=length(k); i++) {
                    print k[i];
                }
            }
        }
    ' ../ES0304_run_finalEMresults.chr"$1".clusters.txt > Clusters_ES0"$2"_as_seed/clusters_final_ES0"$2"_as_seed_chr"$1".txt
}

# Loop through chromosomes
for chr in {1..8}; do
    process_clusters "$chr" 3
    process_clusters "$chr" 4
done
