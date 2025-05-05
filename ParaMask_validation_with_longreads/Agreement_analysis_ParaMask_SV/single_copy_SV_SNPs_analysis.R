#### Script to analyse 'false' single-copy called SNPs by PM overlapping with SVs #####
### 10.04.2025 ###
# define chromosome vector:
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

#### ========= Read in files: ======== ####
#### 1) Read EM calls of ParaMask
for (chr in chromosomes) {
  PM_currecnt_chr_df <- read.delim(paste("ES0304_run_finalEMresults.", chr,".finalClass.het", sep = ""))
  EM_ES03_ES04_finalClass_original[[chr]] <- PM_currecnt_chr_df
}

#### 2) Cluster files of PM
EM_ES03_ES04_cluster_original <- list()
for (chr in chromosomes) {
  PM_currecnt_chr_df <- read.delim(paste("ES0304_run_finalEMresults.", chr,".clusters.txt", sep = ""))
  EM_ES03_ES04_cluster_original[[chr]] <- PM_currecnt_chr_df
}


#### 3) read SV long-read calls
ES03_014_SV_loop_merged <- read.delim("ES03_014_only_overlapping_dup_with_header.txt")
ES04_014_SV_loop_merged <- read.delim("ES04_014_only_overlapping_dup_with_header.txt")



#### ======== filter EM calls of PM for only overlapping SV region =========
# for ES03-014:
PM_calls_subset_all_SV_ES03_list <- list()
for (chr in chromosomes) {
  current_PM_FP  <- EM_ES03_ES04_finalClass_original[[chr]]
  currect_SV_03 <- ES03_014_SV_loop_merged[ES03_014_SV_loop_merged$Chromosome == chr, ]
  current_PM_FP_sub_list <- list()
  for (row_num in 1:nrow(currect_SV_03)) {
    bin_stop <- currect_SV_03$STOP[row_num]
    bin_start <- currect_SV_03$START[row_num]

    current_PM_FP_sub <- current_PM_FP[current_PM_FP$Position >= bin_start & current_PM_FP$Position <= bin_stop,]
    current_PM_FP_sub$SV_number <- rep(row_num, length(current_PM_FP_sub$Chromosome))
    current_PM_FP_sub_list[[row_num]] <- current_PM_FP_sub
  }
  current_PM_FP_data_all_SVs <- bind_rows(current_PM_FP_sub_list)
  PM_calls_subset_all_SV_ES03_list[[chr]] <- current_PM_FP_data_all_SVs
}
length(PM_calls_subset_all_SV_ES03_df$Chromosome)
PM_calls_subset_all_SV_ES03_df$finalClass_label[PM_calls_subset_all_SV_ES03_df$finalClass==1] <- "multicopy"
PM_calls_subset_all_SV_ES03_df$finalClass_label[PM_calls_subset_all_SV_ES03_df$finalClass==0] <- "singlecopy"


# for ES04-014:
PM_calls_subset_all_SV_ES04_list <- list()
for (chr in chromosomes) {
  current_PM_FP  <- EM_ES03_ES04_finalClass_original[[chr]]
  currect_SV_03 <- ES04_014_SV_loop_merged[ES04_014_SV_loop_merged$Chromosome == chr, ]
  current_PM_FP_sub_list <- list()
  for (row_num in 1:nrow(currect_SV_03)) {
    bin_stop <- currect_SV_03$STOP[row_num]
    bin_start <- currect_SV_03$START[row_num]

    current_PM_FP_sub <- current_PM_FP[current_PM_FP$Position >= bin_start & current_PM_FP$Position <= bin_stop,]
    current_PM_FP_sub$SV_number <- rep(row_num, length(current_PM_FP_sub$Chromosome))
    current_PM_FP_sub_list[[row_num]] <- current_PM_FP_sub
  }
  current_PM_FP_data_all_SVs <- bind_rows(current_PM_FP_sub_list)
  PM_calls_subset_all_SV_ES04_list[[chr]] <- current_PM_FP_data_all_SVs
}
# Combine all chromosomes into one final dataframe
PM_calls_subset_all_SV_ES04_df <- dplyr::bind_rows(PM_calls_subset_all_SV_ES04_list)
PM_calls_subset_all_SV_ES04_df$finalClass_label[PM_calls_subset_all_SV_ES04_df$finalClass==1] <- "multicopy"
PM_calls_subset_all_SV_ES04_df$finalClass_label[PM_calls_subset_all_SV_ES04_df$finalClass==0] <- "singlecopy"


#### ======= Analysis 1:  number of het Genotypes for finalClass == 0 of last seed ======== ####
#### Check the numbers of hetGenotypes in the last seed ####
# for ES03
PM_calls_subset_all_SV_ES03_list_cluster_info <- list()

for (chr in chromosomes) {
  cluster_current <- EM_ES03_ES04_cluster_original[[chr]]
  current_PM_FP  <-  PM_calls_subset_all_SV_ES03_list[[chr]]

  # Initialize columns once outside the loop
  current_PM_FP$clusterCause <- NA
  current_PM_FP$CovHetOfLastSeed <- NA
  current_PM_FP$hetGenOfLastSeed <- NA

  for (pos in unique(current_PM_FP$Position)) {
    # Find matching index in cluster_current
    match_idx <- which(cluster_current$position == pos)

    if (length(match_idx) > 0) {
      assign_idx <- which(current_PM_FP$Position == pos)
      current_PM_FP$clusterCause[assign_idx] <- cluster_current$clusterCause[match_idx]
      current_PM_FP$CovHetOfLastSeed[assign_idx] <- cluster_current$CovHetOfLastSeed[match_idx]
      current_PM_FP$hetGenOfLastSeed[assign_idx] <- cluster_current$hetGenOfLastSeed[match_idx]
    }
  }

  PM_calls_subset_all_SV_ES03_list_cluster_info[[chr]] <- current_PM_FP
}

df <- bind_rows(PM_calls_subset_all_SV_ES03_list_cluster_info)
df$hetGenOfLastSeed <- as.character(df$hetGenOfLastSeed)

# Create a new column to hold the parsed count
df$hetGenCountFromLastSeed <- NA_integer_

# Store the hetGenOfLastSeed string of the last seen seed
last_het_gen_string <- NA

for (i in seq_len(nrow(df))) {
  # Check if it's a seed SNP
  if (df$finalClass[i] == 1 && df$clusterCause[i] == "seed") {
    last_het_gen_string <- df$hetGenOfLastSeed[i]

  } else if (df$finalClass[i] == 0) {
    # If it's a finalClass 0, parse the last seed’s hetGen string
    if (!is.na(last_het_gen_string) && nchar(last_het_gen_string) > 0) {
      het_count <- length(unlist(strsplit(last_het_gen_string, ":")))
      df$hetGenCountFromLastSeed[i] <- het_count
    } else {
      df$hetGenCountFromLastSeed[i] <- 0
    }
  }
}
hist(df$hetGenCountFromLastSeed[df$finalClass==0], breaks = 86)



PM_calls_subset_all_SV_ES03_list_cluster_df$hetGen_count <- sapply(
  strsplit(PM_calls_subset_all_SV_ES03_list_cluster_df$hetGenOfLastSeed, ":"),
  length
)
write.table(df, "ES03_FN_SNPs_hetGen_of_last_seed_SNPs.txt")
#df <- read.table("FN_SNPs_hetGen_of_last_seed_SNPs.txt")


# for ES04:
PM_calls_subset_all_SV_ES04_list_cluster_info <- list()

for (chr in chromosomes) {
  cluster_current <- EM_ES03_ES04_cluster_original[[chr]]
  current_PM_FP  <-  PM_calls_subset_all_SV_ES04_list[[chr]]

  # Initialize columns once outside the loop
  current_PM_FP$clusterCause <- NA
  current_PM_FP$CovHetOfLastSeed <- NA
  current_PM_FP$hetGenOfLastSeed <- NA

  for (pos in unique(current_PM_FP$Position)) {
    # Find matching index in cluster_current
    match_idx <- which(cluster_current$position == pos)

    if (length(match_idx) > 0) {
      assign_idx <- which(current_PM_FP$Position == pos)
      current_PM_FP$clusterCause[assign_idx] <- cluster_current$clusterCause[match_idx]
      current_PM_FP$CovHetOfLastSeed[assign_idx] <- cluster_current$CovHetOfLastSeed[match_idx]
      current_PM_FP$hetGenOfLastSeed[assign_idx] <- cluster_current$hetGenOfLastSeed[match_idx]
    }
  }

  PM_calls_subset_all_SV_ES04_list_cluster_info[[chr]] <- current_PM_FP
}

df_ES04 <- bind_rows(PM_calls_subset_all_SV_ES04_list_cluster_info)
df_ES04$hetGenOfLastSeed <- as.character(df_ES04$hetGenOfLastSeed)

# Create a new column to hold the parsed count
df_ES04$hetGenCountFromLastSeed <- NA_integer_
last_het_gen_string <- NA

for (i in seq_len(nrow(df_ES04))) {
  # Check if it's a seed SNP
  if (df_ES04$finalClass[i] == 1 && df_ES04$clusterCause[i] == "seed") {
    last_het_gen_string <- df_ES04$hetGenOfLastSeed[i]

  } else if (df_ES04$finalClass[i] == 0) {
    # If it's a finalClass 0, parse the last seed’s hetGen string
    if (!is.na(last_het_gen_string) && nchar(last_het_gen_string) > 0) {
      het_count <- length(unlist(strsplit(last_het_gen_string, ":")))
      df_ES04$hetGenCountFromLastSeed[i] <- het_count
    } else {
      df_ES04$hetGenCountFromLastSeed[i] <- 0
    }
  }
}


write.table(df_ES04, "ES04_FN_SNPs_hetGen_of_last_seed_SNPs.txt")



#### ======== Analysis 2: Check percentage of SNPs across SV windows 0f 0.1 length SV ========= #####

FP_SNP_analysis_list_ES03 <- list()

## for ES03-014:
for (chr in chromosomes) {
  current_PM_FP  <- EM_ES03_ES04_finalClass_original[[chr]][EM_ES03_ES04_finalClass_original[[chr]]$finalClass==0,]
  currect_SV_03 <- ES03_014_SV_loop_merged[ES03_014_SV_loop_merged$Chromosome == chr, ]

  # Reinitialize for each chromosome
  tenpercent_SV <- data.frame(
    Chromosome = character(0),
    SV_number = integer(0),
    START = integer(0),
    STOP = integer(0),
    Percentage_Bin = character(0),
    Number_of_FP_SNPs = integer(0),
    Proportion_of_FP_SNP_across_all_FP = numeric(0)
  )

  for (row_num in 1:nrow(currect_SV_03)) {
    ten_percent <- (currect_SV_03$STOP[row_num] - currect_SV_03$START[row_num]) * 0.1

    for (i in 1:10) {
      bin_start <- currect_SV_03$START[row_num] + (i - 1) * ten_percent
      bin_stop <- bin_start + ten_percent
      percentage_label <- paste0((i - 1) * 10, "_", i * 10)

      tenpercent_SV <- rbind(tenpercent_SV, data.frame(
        Chromosome = chr,
        SV_number = row_num,
        START = round(bin_start),
        STOP = round(bin_stop),
        Percentage_Bin = percentage_label,
        Number_of_FP_SNPs = NA,
        Proportion_of_FP_SNP_across_all_FP = NA
      ))
    }
  }

  # Compute number and proportion of FP SNPs per bin
  for (bin_row in 1:nrow(tenpercent_SV)) {
    bin_start <- tenpercent_SV$START[bin_row]
    bin_stop <- tenpercent_SV$STOP[bin_row]
    row_num <- tenpercent_SV$SV_number[bin_row]

    # Count FP SNPs in current bin
    tenpercent_SV$Number_of_FP_SNPs[bin_row] <- sum(
      current_PM_FP$Position >= bin_start & current_PM_FP$Position <= bin_stop
    )

    # Count all FP SNPs in full SV
    all_FP_SNPs <- sum(
      current_PM_FP$Position >= currect_SV_03$START[row_num] &
        current_PM_FP$Position <= currect_SV_03$STOP[row_num]
    )

    if (all_FP_SNPs > 0) {
      tenpercent_SV$Proportion_of_FP_SNP_across_all_FP[bin_row] <-
        tenpercent_SV$Number_of_FP_SNPs[bin_row] / all_FP_SNPs
    } else {
      tenpercent_SV$Proportion_of_FP_SNP_across_all_FP[bin_row] <- NA
    }
  }

  FP_SNP_analysis_list_ES03[[chr]] <- tenpercent_SV
}

# Combine all chromosomes into one final dataframe
FP_SNP_analysis_df_ES03 <- dplyr::bind_rows(FP_SNP_analysis_list_ES03)



#### ES04:
FP_SNP_analysis_list_ES04 <- list()

for (chr in chromosomes) {
  current_PM_FP  <- EM_ES03_ES04_finalClass_original[[chr]][EM_ES03_ES04_finalClass_original[[chr]]$finalClass==0,]
  currect_SV_04 <- ES04_014_SV_loop_merged[ES04_014_SV_loop_merged$Chromosome == chr, ]

  # Reinitialize for each chromosome
  tenpercent_SV <- data.frame(
    Chromosome = character(0),
    SV_number = integer(0),
    START = integer(0),
    STOP = integer(0),
    Percentage_Bin = character(0),
    Number_of_FP_SNPs = integer(0),
    Proportion_of_FP_SNP_across_all_FP = numeric(0)
  )

  for (row_num in 1:nrow(currect_SV_04)) {
    ten_percent <- (currect_SV_04$STOP[row_num] - currect_SV_04$START[row_num]) * 0.1

    for (i in 1:10) {
      bin_start <- currect_SV_04$START[row_num] + (i - 1) * ten_percent
      bin_stop <- bin_start + ten_percent
      percentage_label <- paste0((i - 1) * 10, "_", i * 10)

      tenpercent_SV <- rbind(tenpercent_SV, data.frame(
        Chromosome = chr,
        SV_number = row_num,
        START = round(bin_start),
        STOP = round(bin_stop),
        Percentage_Bin = percentage_label,
        Number_of_FP_SNPs = NA,
        Proportion_of_FP_SNP_across_all_FP = NA
      ))
    }
  }

  # Compute number and proportion of FP SNPs per bin
  for (bin_row in 1:nrow(tenpercent_SV)) {
    bin_start <- tenpercent_SV$START[bin_row]
    bin_stop <- tenpercent_SV$STOP[bin_row]
    row_num <- tenpercent_SV$SV_number[bin_row]

    # Count FP SNPs in current bin
    tenpercent_SV$Number_of_FP_SNPs[bin_row] <- sum(
      current_PM_FP$Position >= bin_start & current_PM_FP$Position <= bin_stop
    )

    # Count all FP SNPs in full SV
    all_FP_SNPs <- sum(
      current_PM_FP$Position >= currect_SV_04$START[row_num] &
        current_PM_FP$Position <= currect_SV_04$STOP[row_num]
    )

    if (all_FP_SNPs > 0) {
      tenpercent_SV$Proportion_of_FP_SNP_across_all_FP[bin_row] <-
        tenpercent_SV$Number_of_FP_SNPs[bin_row] / all_FP_SNPs
    } else {
      tenpercent_SV$Proportion_of_FP_SNP_across_all_FP[bin_row] <- NA
    }
  }

  FP_SNP_analysis_list_ES04[[chr]] <- tenpercent_SV
}

# Combine all chromosomes into one final dataframe
FP_SNP_analysis_df_ES04 <- dplyr::bind_rows(FP_SNP_analysis_list_ES04)
sum(FP_SNP_analysis_df_ES04$Number_of_FP_SNPs)
sum(FP_SNP_analysis_df_ES03$Number_of_FP_SNPs)

#######

FP_plot_df_03 <- FP_SNP_analysis_df_ES03[!is.na(FP_SNP_analysis_df_ES03$Proportion_of_FP_SNP_across_all_FP), ]
FP_plot_df_04 <- FP_SNP_analysis_df_ES04[!is.na(FP_SNP_analysis_df_ES04$Proportion_of_FP_SNP_across_all_FP), ]

sum(FP_plot_df_03$Number_of_FP_SNPs)
sum(FP_plot_df_04$Number_of_FP_SNPs)

sum(FP_plot_df_03$Proportion_of_FP_SNP_across_all_FP)
sum(FP_plot_df_04$Proportion_of_FP_SNP_across_all_FP)

for (chr in chromosomes) {
  chr_sub <- FP_plot_df_04[FP_plot_df_04$Chromosome==chr,]
  for (SV in unique(chr_sub$SV_number)) {
    sub <- chr_sub[chr_sub$SV_number==SV,]
    print(sum(sub$Proportion_of_FP_SNP_across_all_FP))
  }
}


ES03_false_pos <- ggplot(FP_plot_df_03, aes(x = Percentage_Bin, y = Proportion_of_FP_SNP_across_all_FP)) +
  geom_quasirandom(shape = 20, color = "black", size = 4) +
  geom_boxplot(fill = NA, size = 1.2) +
  labs(
    title = "ES03: 'False' single-copy SNP distribution across duplications (10% windows) ",
    x = "SV Bin (% along SV)",
    y = "Proportion of FP SNPs (across all FP SNPs in SV)",
  ) +
  theme(
    text = element_text(size = 13),
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15, vjust = 0.1),
    axis.title.y = element_text(size = 15),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(fill = NA, color = "white"),
    legend.position = "none"
  )


ES04_false_pos <- ggplot(FP_plot_df_04, aes(x = Percentage_Bin, y = Proportion_of_FP_SNP_across_all_FP)) +
  geom_quasirandom(shape = 20, color = "black", size = 4) +
  geom_boxplot(fill = NA, size = 1.2) +
  labs(
    title = "ES04: 'False' single-copy SNP distribution across duplications (10% windows) ",
    x = "SV Bin (% along SV)",
    y = "Proportion of FP SNPs (across all FP SNPs in SV)",
  ) +
  theme(
    text = element_text(size = 13),
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 15, vjust = 0.1),
    axis.title.y = element_text(size = 15),
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(fill = NA, color = "white"),
    legend.position = "none"
  )


ggarrange(ES03_false_pos,
          ES04_false_pos,
          nrow=2)
