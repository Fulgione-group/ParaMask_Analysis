# Load required libraries
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggbeeswarm)

# Define paths
base_path <- "~/PATH/"
vcf_path <- file.path(base_path, "/Path/to/SV/")
final_path <- file.path(base_path, "Path/to/PM_calls")

# Read in the ES03_014_SV_loop_merged file
ES03_014_SV_loop_merged <- read.delim(file.path(vcf_path, "ES03_014_SV_file.txt"))

# Define chromosomes
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

# Initialize lists
rCNV_paralog_list_original <- list()
EM_ES03_ES04_finalClass_original <- list()
EM_ES03_ES04_finalClass_filtered <- list()
rCNV_paralog_list_filtered <- list()
Merged_list_rCNV_PM <- list()

# Loop through chromosomes and process data
for (chr in chromosomes) {

  # Read all ParaMask results
  PM_currecnt_chr_df <- read.delim(file.path(final_path, paste("ES0304_run_finalEMresults.", chr, ".finalClass.het", sep = "")))
  EM_ES03_ES04_finalClass_original[[chr]] <- PM_currecnt_chr_df

  # Read all rCNV results
  rCNV_currecnt_chr_df <- read.csv(file.path("Y:/dep_coupland/grp_fulgione/helene/data/rCNV", chr, "cnv.txt"), sep = "")
  colnames(rCNV_currecnt_chr_df) <- c("chromosome", "Position", "ID", "NHet", "medRatio", "NHomRef", "NHomAlt", "propHomAlt", "Nsamp", "dup.stat", "class")

  # Assign finalClass column
  rCNV_currecnt_chr_df$finalClass <- rep(0, length(rCNV_currecnt_chr_df$chromosome))

  # Assign final class 1 or 0
  rCNV_currecnt_chr_df$finalClass[rCNV_currecnt_chr_df$class == "cnv"] <- 1

  rCNV_paralog_list_original[[chr]] <- rCNV_currecnt_chr_df

  merged_df <- merge(rCNV_currecnt_chr_df, PM_currecnt_chr_df, by = c("Position"), all = TRUE, suffixes = c("_rCNV", "_EM"))
  # Replace NA values in finalClass_rCNV with 0
  merged_df$finalClass_rCNV[is.na(merged_df$finalClass_rCNV)] <- 0
  merged_df$class[is.na(merged_df$class)] <- "uncertain"
  merged_df$chromosome[is.na(merged_df$chromosome)] <- chr

  # Rename columns to match original rCNV_original_table structure, if necessary
  colnames(merged_df)[colnames(merged_df) == "finalClass_rCNV"] <- "finalClass"
  colnames(merged_df)[colnames(merged_df) == "finalClass_EM"] <- "finalClass_EM"

  # Optionally, drop columns from EM_ES03_ES04_finalClass_original_table that are not needed
  # This step depends on your specific requirements
  columns_to_keep <- colnames(rCNV_currecnt_chr_df)
  merged_df_update <- merged_df[, columns_to_keep]

  Merged_list_rCNV_PM[[chr]] <- merged_df_update
  # Find overlapping positions
  overlapping_positions <- intersect(PM_currecnt_chr_df$Position, rCNV_currecnt_chr_df$Position)

  # Filter the data frames to include only overlapping positions
  filtered_em_df <- PM_currecnt_chr_df[PM_currecnt_chr_df$Position %in% overlapping_positions, ]
  filtered_rCNV_df <- rCNV_currecnt_chr_df[rCNV_currecnt_chr_df$Position %in% overlapping_positions, ]

  EM_ES03_ES04_finalClass_filtered[[chr]] <- filtered_em_df
  rCNV_paralog_list_filtered[[chr]] <- filtered_rCNV_df
}
Merged_table_rCNV_PM <- bind_rows(Merged_list_rCNV_PM)
table(Merged_table_rCNV_PM$class)

## Original Table ####
# Different positions
chr <- "chr1"
unique_rCNV_positions <- setdiff(EM_ES03_ES04_finalClass_original[[chr]]$Position, rCNV_paralog_list_original[[chr]]$Position)
unique_rCNV_positions

## STATS PM
EM_ES03_ES04_finalClass_original_table <- bind_rows(EM_ES03_ES04_finalClass_original)
table(EM_ES03_ES04_finalClass_original_table$finalClass)
EM_ES03_ES04_finalClass_table_filtered <- bind_rows(EM_ES03_ES04_finalClass_filtered)
length(EM_ES03_ES04_finalClass_table_filtered$finalClass)

## FILTERED POSITION TABLE ####
## OVERLAP METHODS ####
# uncertain SNPs included
Agreement_table_rCNV_list_uncertain_included <- list()
for (chr in chromosomes) {

  # Compare the finalClass values
  Agreement_table_rCNV_uncertain_included <- data.frame(
    Chromosome <- EM_ES03_ES04_finalClass_original[[chr]]$Chromosome,
    Position <- EM_ES03_ES04_finalClass_original[[chr]]$Position,
    MAF <- EM_ES03_ES04_finalClass_original[[chr]]$Minor.allele.freq,
    ParaMask_fC <- EM_ES03_ES04_finalClass_original[[chr]]$finalClass,
    rCNV_fC <- Merged_list_rCNV_PM[[chr]]$finalClass,
    rCNV_class <- Merged_list_rCNV_PM[[chr]]$class
  )
  colnames(Agreement_table_rCNV_uncertain_included) <- c("Chromosome", "Position","MAF", "ParaMask_fC", "rCNV_fC", "rCNV_class")

  Agreement_table_rCNV_uncertain_included$Overlap_PM_rCNV <- rep(FALSE, length(Agreement_table_rCNV_uncertain_included$Chromosome))

  for (row_num in 1:length(Agreement_table_rCNV_uncertain_included$Chromosome)) {
    if (Agreement_table_rCNV_uncertain_included$ParaMask_fC[row_num] == Agreement_table_rCNV_uncertain_included$rCNV_fC[row_num]) {
      Agreement_table_rCNV_uncertain_included$Overlap_PM_rCNV[row_num] <- TRUE
    }
  }
  Agreement_table_rCNV_list_uncertain_included[[chr]] <- Agreement_table_rCNV_uncertain_included
}
Agreement_table_rCNV_all_chr_uncertain_included <- bind_rows(Agreement_table_rCNV_list_uncertain_included)
Agreement_table_rCNV_all_chr_uncertain_included$Overlap_PM_rCNV[Agreement_table_rCNV_all_chr_uncertain_included$rCNV_class == "uncertain"] <- "uncertain"

########################################
# Only overlapping positions
Agreement_table_rCNV_list_only_overlapping <- list()
for (chr in chromosomes) {

  # Compare the finalClass values
  Agreement_table_rCNV_only_overlapping <- data.frame(
    Chromosome <- EM_ES03_ES04_finalClass_filtered[[chr]]$Chromosome,
    Position <- EM_ES03_ES04_finalClass_filtered[[chr]]$Position,
    MAF <- EM_ES03_ES04_finalClass_filtered[[chr]]$Minor.allele.freq,
    ParaMask_fC <- EM_ES03_ES04_finalClass_filtered[[chr]]$finalClass,
    rCNV_fC <- rCNV_paralog_list_filtered[[chr]]$finalClass,
    rCNV_class <- rCNV_paralog_list_filtered[[chr]]$class
  )
  colnames(Agreement_table_rCNV_only_overlapping) <- c("Chromosome", "Position","MAF", "ParaMask_fC", "rCNV_fC", "rCNV_class")

  Agreement_table_rCNV_only_overlapping$Overlap_PM_rCNV <- rep(FALSE, length(Agreement_table_rCNV_only_overlapping$Chromosome))

  for (row_num in 1:length(Agreement_table_rCNV_only_overlapping$Chromosome)) {
    if (Agreement_table_rCNV_only_overlapping$ParaMask_fC[row_num] == Agreement_table_rCNV_only_overlapping$rCNV_fC[row_num]) {
      Agreement_table_rCNV_only_overlapping$Overlap_PM_rCNV[row_num] <- TRUE
    }
  }
  Agreement_table_rCNV_list_only_overlapping[[chr]] <- Agreement_table_rCNV_only_overlapping
}
Agreement_table_rCNV_all_chr_only_overlapping <- bind_rows(Agreement_table_rCNV_list_only_overlapping)

length(Agreement_table_rCNV_all_chr_only_overlapping$Chromosome)

# STATS FILTERED POSITIONS ####
length(Agreement_table_rCNV_all_chr_uncertain_included$Chromosome)
length(Agreement_table_rCNV_all_chr_uncertain_included$Chromosome[Agreement_table_rCNV_all_chr_uncertain_included$Overlap_PM_rCNV == "TRUE"])
table(Agreement_table_rCNV_all_chr_uncertain_included$Overlap_PM_rCNV)
table(Agreement_table_rCNV_all_chr_uncertain_included$Overlap_PM_rCNV[Agreement_table_rCNV_all_chr_uncertain_included$rCNV_class != "uncertain"])

# OVERLAP with SV
ES03_014_SV_loop_merged <- read.delim(file.path(vcf_path, "ES03_014_only_overlapping_dup_with_header.txt"))
ES04_014_SV_loop_merged <- read.delim(file.path(vcf_path, "ES04_014_only_overlapping_dup_with_header.txt"))

rCNV_finalClass_SV_merged_list <- list()
for (chr in chromosomes) {

  current_chr_df_rCNV <- Agreement_table_rCNV_list_only_overlapping[[chr]]
  # Replacing columns that are not needed from SV table
  SV_col_names <- colnames(current_chr_df_rCNV)

  currect_SV_03 <- ES03_014_SV_loop_merged[, c(1, 2, 3, 4, 10)]
  currect_SV_03 <- currect_SV_03[currect_SV_03$CHR == chr, ]
  currect_SV_04 <- ES04_014_SV_loop_merged[, c(1, 2, 3, 4, 10)]
  currect_SV_04 <- currect_SV_04[currect_SV_04$CHR == chr, ]

  # Initiating columns
  current_chr_df_rCNV$overlap_SV_03 <- rep(FALSE, length(current_chr_df_rCNV$Chromosome))
  current_chr_df_rCNV$overlap_SV_04 <- rep(FALSE, length(current_chr_df_rCNV$Chromosome))
  current_chr_df_rCNV$mean_fC_SV_03 <- rep(NA, length(current_chr_df_rCNV$Chromosome))
  current_chr_df_rCNV$mean_fC_SV_04 <- rep(NA, length(current_chr_df_rCNV$Chromosome))

  for (row_num in 1:length(current_chr_df_rCNV$Chromosome)) {

    SV_03_indices <- which(currect_SV_03$start.pos <= current_chr_df_rCNV$Position[row_num] & currect_SV_03$end.pos >= current_chr_df_rCNV$Position[row_num])
    SV_04_indices <- which(currect_SV_04$start.pos <= current_chr_df_rCNV$Position[row_num] & currect_SV_04$end.pos >= current_chr_df_rCNV$Position[row_num])

    if (length(SV_03_indices) > 0) {
      current_chr_df_rCNV$overlap_SV_03[row_num] <- TRUE
      current_chr_df_rCNV$mean_fC_SV_03[row_num] <- mean(currect_SV_03$log2.copyRatio[SV_03_indices])
    }

    if (length(SV_04_indices) > 0) {
      current_chr_df_rCNV$overlap_SV_04[row_num] <- TRUE
      current_chr_df_rCNV$mean_fC_SV_04[row_num] <- mean(currect_SV_04$log2.copyRatio[SV_04_indices])
    }
  }

  rCNV_finalClass_SV_merged_list[[chr]] <- current_chr_df_rCNV
}
rCNV_finalClass_SV_merged_table <- bind_rows(rCNV_finalClass_SV_merged_list)
