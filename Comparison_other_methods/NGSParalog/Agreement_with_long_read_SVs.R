### OVERLAP WITH LONG READS ###


ES03_014_SV_loop_merged <- read.delim("~/PATH/TO/SV/ES03_014_only_overlapping_dup_with_header.txt")

ES04_014_SV_loop_merged <- read.delim("~/PATH/TO/SV/ES04_014_only_overlapping_dup_with_header.txt")


ES04_014_overlapping_dup_list <- list()
ES03_014_overlapping_dup_list <- list()
Agreement_list_with_SVs <- list()

## CHeck Overlap with SVs ##
for (chr in chromosomes) {
  #Grep dataframe of corresponding chromosome and create new columns 
  Agreement_current_chr <- Agreement_table_all_chr[Agreement_table_all_chr$Chromosome==chr,]
  
  SV_current_chr_ES03_014 <- ES03_014_SV_loop_merged[ES03_014_SV_loop_merged$Chromosome==chr,]
  
  SV_current_chr_ES03_014$mean_final_class_ngs_BH <- rep(0, length(SV_current_chr_ES03_014$Chromosome))
  SV_current_chr_ES03_014$mean_final_class_ngs_bonferroni <- rep(0, length(SV_current_chr_ES03_014$Chromosome))
  
  SV_current_chr_ES04_014 <- ES04_014_SV_loop_merged[ES04_014_SV_loop_merged$Chromosome==chr,]
  
  SV_current_chr_ES04_014$mean_final_class_ngs_BH <- rep(0, length(SV_current_chr_ES04_014$Chromosome))
  SV_current_chr_ES04_014$mean_final_class_ngs_bonferroni <- rep(0, length(SV_current_chr_ES04_014$Chromosome))
  
  Agreement_current_chr$SV_call_ES03_014 <- rep("no", length(Agreement_current_chr$Chromosome))
  Agreement_current_chr$SV_call_fC_ES03_014 <- rep(0, length(Agreement_current_chr$Chromosome))
  
  Agreement_current_chr$SV_call_ES04_014 <- rep("no", length(Agreement_current_chr$Chromosome))
  Agreement_current_chr$SV_call_fC_ES04_014 <- rep(0, length(Agreement_current_chr$Chromosome))
  
  # Assign SV information of ES03-014
  for (row_num in 1:length(SV_current_chr_ES03_014$Chromosome)) {
    Agreement_current_chr$SV_call_fC_ES03_014[Agreement_current_chr$Position >= SV_current_chr_ES03_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES03_014$STOP[row_num]] <- 1
    Agreement_current_chr$SV_call_ES03_014[Agreement_current_chr$Position >= SV_current_chr_ES03_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES03_014$STOP[row_num]] <- "yes"
    print("ES03_014")
    print(mean(Agreement_current_chr$NGSparalog_fC[Agreement_current_chr$Position >= SV_current_chr_ES03_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES03_014$STOP[row_num]]))
    # Benjamin Hochberg
    SV_current_chr_ES03_014$mean_final_class_ngs_BH[row_num] <- mean(Agreement_current_chr$NGSparalog_fC_BH[Agreement_current_chr$Position >= SV_current_chr_ES03_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES03_014$STOP[row_num]])
    # Bonferroni
    SV_current_chr_ES03_014$mean_final_class_ngs_bonferroni[row_num] <- mean(Agreement_current_chr$NGSparalog_fC_bonferroni[Agreement_current_chr$Position >= SV_current_chr_ES03_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES03_014$STOP[row_num]])
  }
  
  # Assign SV information of ES04-014
  for (row_num in 1:length(SV_current_chr_ES04_014$Chromosome)) {
    Agreement_current_chr$SV_call_fC_ES04_014[Agreement_current_chr$Position >= SV_current_chr_ES04_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES04_014$STOP[row_num]] <- 1
    Agreement_current_chr$SV_call_ES04_014[Agreement_current_chr$Position >= SV_current_chr_ES04_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES04_014$STOP[row_num]] <- "yes"
    print("ES04_014")
    print( mean(Agreement_current_chr$NGSparalog_fC[Agreement_current_chr$Position >= SV_current_chr_ES04_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES04_014$STOP[row_num]]))
    
    SV_current_chr_ES04_014$mean_final_class_ngs_BH[row_num] <- mean(Agreement_current_chr$NGSparalog_fC_BH[Agreement_current_chr$Position >= SV_current_chr_ES04_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES04_014$STOP[row_num]])
    SV_current_chr_ES04_014$mean_final_class_ngs_bonferroni[row_num] <- mean(Agreement_current_chr$NGSparalog_fC_bonferroni[Agreement_current_chr$Position >= SV_current_chr_ES04_014$START[row_num] &  Agreement_current_chr$Position <= SV_current_chr_ES04_014$STOP[row_num]])
  }
  
  ES03_014_overlapping_dup_list[[chr]] <- SV_current_chr_ES03_014
  ES04_014_overlapping_dup_list[[chr]] <- SV_current_chr_ES04_014
  Agreement_list_with_SVs[[chr]] <- Agreement_current_chr
}
