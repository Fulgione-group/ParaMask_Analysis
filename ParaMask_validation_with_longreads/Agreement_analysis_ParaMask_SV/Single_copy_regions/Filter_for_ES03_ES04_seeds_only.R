## RStudio: Read in files (.cluster.txt; .het; .bed)

## READ IN ParaMask output (chr1)
# Store file paths in a list
file_paths_het <- paste0("~/PATH/TO/PARAMASK/HET-FILES/ES0304_run_finalEMresults.chr", chromosomes, ".finalClass.het")

# Read HET files into a list
het_files <- lapply(file_paths_het, read.delim)

# Store cluster file paths in a list
cluster_file_paths <- paste0("~/PATH/TO/PARAMASK/CLUSTER-FILTERED-FILES/clusters_final_ES03_as_seed_chr", chromosomes, ".txt")

# Read cluster files into a list
cluster_files <- lapply(cluster_file_paths, read.delim, header = FALSE)

# Initialize empty lists
Clusters_final_ES03_as_seed_list <- list()
ES03_014_EM_finalClass_filtered <- list()

# Loop through the chromosomes
for (i in seq_along(chromosomes)) {
  # Extract data frames from the lists
  subset_het_filtered <- het_files[[i]]
  subset_cluster_filtered <- cluster_files[[i]]

  # Read bed file
  current_bed_file <- read.delim(paste0("~/PATH/TO/PARAMASK/BED-FILES/ES0304_run_finalEMresults.", chromosomes[i], ".finalClass.bed"), header = TRUE)

  # Extract rows with paralogs
  current_paralog_df <- current_bed_file[current_bed_file$type.0.single.copy.1.multi.copy == 1, ]

  # Filter for size
  size_filter_df <- current_paralog_df[current_paralog_df$End - current_paralog_df$Start < 35, ]

  # Filter clusters based on size
  subset_cluster_filtered_size <- subset_cluster_filtered[!(subset_cluster_filtered$V3 %in% size_filter_df$Cluster), ]

  # Update finalClass in subset_het_filtered
  subset_het_filtered$finalClass[subset_het_filtered$Position %in% subset_cluster_filtered_size$V2] <- 1

  # Add to lists
  ES03_014_EM_finalClass_filtered[[chromosomes[i]]] <- subset_het_filtered
  Clusters_final_ES03_as_seed_list[[chromosomes[i]]] <- subset_cluster_filtered_size
}

# Bind rows of ES03_014_EM_finalClass_filtered
ES03_014_EM_finalClass_filtered_data_frame <- do.call(rbind, ES03_014_EM_finalClass_filtered)

# Count finalClass == 1
sum(ES03_014_EM_finalClass_filtered_data_frame$finalClass == 1)
