## Script to compare overlap between ngsparalog and simulations: ####
### 01.08.2025 ###
## define vecctor with conditions and parameters:
coverage_samples <- "10X_100N"        # Adjust for 5X_15N
condition_vector <- c("HW_Rep1", "HW_Rep2", "HW_Rep3", "Fis0.9_Rep1", "Fis0.9_Rep2", "Fis0.9_Rep3")

## list to store results for each condition ##
results_bonferroni <- list()
for (condition in condition_vector) {
  # Load files
  print(condition)
  position_file <- read.delim(paste0("N:/path/to/Paramask_sim_fastas/Position_files_", condition, ".txt"), header=FALSE)
  ngsparalog_output <- read.delim(paste0("N:/path/to/ParaMask/ngsParalog/seq_sim/",coverage_samples,"/", condition, "/", condition, "_ngsParalog_output.lr"), header=FALSE)
  # Prepare dataframe
  colnames(ngsparalog_output) <- c("chromosome", "Position", "LR_1", "LR_2", "final_LR")
  ngsparalog_output$pval <- 0.5 * pchisq(ngsparalog_output$final_LR, df=1, lower.tail=FALSE)
  ngsparalog_output$pval.adj_BH <- p.adjust(ngsparalog_output$pval, method="BH")
  ngsparalog_output$pval.adj_bonferroni <- p.adjust(ngsparalog_output$pval, method="bonferroni")

  # Assign class labels
  ngsparalog_output$finalClass_BH <- 0
  ngsparalog_output$finalClass_bonferroni <- 0
  ngsparalog_output$finalClass_BH[ngsparalog_output$pval.adj_BH < 0.001] <- 1
  ngsparalog_output$finalClass_bonferroni[ngsparalog_output$pval.adj_bonferroni < 0.05] <- 1

  # Filter overlaps
  overlapping_positions <- intersect(position_file$V1, ngsparalog_output$Position)
  filtered_ngs_df <- ngsparalog_output[ngsparalog_output$Position %in% overlapping_positions, ]
  filtered_sim_df <- position_file[position_file$V1 %in% overlapping_positions, ]
  print(length(overlapping_positions))
  # Separate by copy type
  multicopy_positions <- position_file$V1[position_file$V2 == "multi-copy"]
  singlecopy_positions <- position_file$V1[position_file$V2 == "single-copy"]

  # Bonferroni overlaps
  mc_df <- filtered_ngs_df[filtered_ngs_df$Position %in% multicopy_positions, ]
  sc_df <- filtered_ngs_df[filtered_ngs_df$Position %in% singlecopy_positions, ]

 # multi_multi <- round(sum(mc_df$finalClass_bonferroni == 1) / length(multicopy_positions), digits = 2)
  multi_multi <- round(sum(mc_df$finalClass_bonferroni == 1) / length(mc_df$finalClass_bonferroni), digits = 2)
  multi_sc <- round(sum(mc_df$finalClass_bonferroni == 0) / length(mc_df$finalClass_bonferroni), digits = 2)
  sc_sc <- round(sum(sc_df$finalClass_bonferroni == 0) / length(sc_df$chromosome), digits = 2)
  sc_mc <- round(sum(sc_df$finalClass_bonferroni == 1) / length(sc_df$chromosome), digits = 2)

  recall <- (
    length(intersect(filtered_sim_df$V1[filtered_sim_df$V2 == "single-copy"], filtered_ngs_df$Position[filtered_ngs_df$finalClass_bonferroni == 0])) +
      length(intersect(filtered_sim_df$V1[filtered_sim_df$V2 == "multi-copy"], filtered_ngs_df$Position[filtered_ngs_df$finalClass_bonferroni == 1]))
  ) / nrow(filtered_ngs_df)
  recall <- round(recall, digits = 2)
  # Store result
  results_bonferroni[[condition]] <- data.frame(
    padjust = "bonferroni",
    dataset = condition,
    coverage_samples = coverage_samples,
    multi_multi = multi_multi,
    multi_sc = multi_sc,
    sc_sc = sc_sc,
    sc_mc = sc_mc,
    recall = recall
  )
}

# Combine all into one dataframe
final_bonferroni_results <- do.call(rbind, results_bonferroni)



### ============= For BH correction ============ ###
results_BH <- list()  # list to collect all results
for (condition in condition_vector) {
  # Load files
  position_file <- read.delim(paste0("N:/path/to/Paramask_sim_fastas/Position_files_", condition, ".txt"), header=FALSE)
  ngsparalog_output <- read.delim(paste0("N:/path/to/ParaMask/ngsParalog/seq_sim/",coverage_samples,"/", condition, "/", condition, ".lr"), header=FALSE)

  # Prepare dataframe
  colnames(ngsparalog_output) <- c("chromosome", "Position", "LR_1", "LR_2", "final_LR")
  ngsparalog_output$pval <- 0.5 * pchisq(ngsparalog_output$final_LR, df=1, lower.tail=FALSE)
  ngsparalog_output$pval.adj_BH <- p.adjust(ngsparalog_output$pval, method="BH")
  ngsparalog_output$pval.adj_bonferroni <- p.adjust(ngsparalog_output$pval, method="bonferroni")

  # Assign class labels
  ngsparalog_output$finalClass_BH <- 0
  ngsparalog_output$finalClass_bonferroni <- 0
  ngsparalog_output$finalClass_BH[ngsparalog_output$pval.adj_BH < 0.001] <- 1
  ngsparalog_output$finalClass_bonferroni[ngsparalog_output$pval.adj_bonferroni < 0.05] <- 1

  # Filter overlaps
  overlapping_positions <- intersect(position_file$V1, ngsparalog_output$Position)
  filtered_ngs_df <- ngsparalog_output[ngsparalog_output$Position %in% overlapping_positions, ]
  filtered_sim_df <- position_file[position_file$V1 %in% overlapping_positions, ]

  # Separate by copy type
  multicopy_positions <- position_file$V1[position_file$V2 == "multi-copy"]
  singlecopy_positions <- position_file$V1[position_file$V2 == "single-copy"]

  # Bonferroni overlaps
  mc_df <- filtered_ngs_df[filtered_ngs_df$Position %in% multicopy_positions, ]
  sc_df <- filtered_ngs_df[filtered_ngs_df$Position %in% singlecopy_positions, ]

  multi_multi <- round(sum(mc_df$finalClass_BH == 1) / length(mc_df$finalClass_BH), digits = 2)
  multi_sc <- round(sum(mc_df$finalClass_BH == 0) / length(mc_df$finalClass_BH), digits = 2)
  sc_sc <- round(sum(sc_df$finalClass_BH == 0) / length(sc_df$chromosome), digits = 2)
  sc_mc <- round(sum(sc_df$finalClass_BH == 1) / length(sc_df$chromosome), digits = 2)

  recall <- (
    length(intersect(filtered_sim_df$V1[filtered_sim_df$V2 == "single-copy"], filtered_ngs_df$Position[filtered_ngs_df$finalClass_BH == 0])) +
      length(intersect(filtered_sim_df$V1[filtered_sim_df$V2 == "multi-copy"], filtered_ngs_df$Position[filtered_ngs_df$finalClass_BH == 1]))
  ) / nrow(filtered_ngs_df)
  recall <- round(recall, digits = 2)
  # Store result
  results_BH[[condition]] <- data.frame(
    padjust = "BH",
    dataset = condition,
    coverage_samples = coverage_samples,
    multi_multi = multi_multi,
    multi_sc = multi_sc,
    sc_sc = sc_sc,
    sc_mc = sc_mc,
    recall = recall
  )
}
# Combine all into one dataframe
final_BH_results <- do.call(rbind, results_BH)
