## OVERLAP METHODS ####

Agreement_table_list <- list()
for (chr in chromosomes) {

    # Compare the finalClass values
  Agreement_table <- data.frame(
    Chromosome <- EM_ES03_ES04_finalClass_filtered[[chr]]$Chromosome,
    Position <- EM_ES03_ES04_finalClass_filtered[[chr]]$Position,
    MAF <- EM_ES03_ES04_finalClass_filtered[[chr]]$Minor.allele.freq,
    ParaMask_fC <- EM_ES03_ES04_finalClass_filtered[[chr]]$finalClass,
    NGSparalog_fC_BH <- NGS_paralog_list_filtered[[chr]]$finalClass_BH,
    NGSparalog_fC_bonferroni <- NGS_paralog_list_filtered[[chr]]$finalClass_bonferroni
    
  )
  colnames(Agreement_table) <- c("Chromosome", "Position","MAF", "ParaMask_fC", "NGSparalog_fC_BH", "NGSparalog_fC_bonferroni")


  Agreement_table$Overlap_BH <- rep(FALSE, length(Agreement_table$Chromosome))
  Agreement_table$Overlap_bonferroni <- rep(FALSE, length(Agreement_table$Chromosome))
  
  for (row_num in 1:length(Agreement_table$Chromosome)) {
    if (Agreement_table$ParaMask_fC[row_num]==Agreement_table$NGSparalog_fC_BH[row_num]) {
      Agreement_table$Overlap_BH[row_num] <- TRUE
    }
    if (Agreement_table$ParaMask_fC[row_num]==Agreement_table$NGSparalog_fC_bonferroni[row_num]) {
      Agreement_table$Overlap_bonferroni[row_num] <- TRUE
    }
    
  }
  Agreement_table_list[[chr]] <- Agreement_table
}
Agreement_table_all_chr <- bind_rows(Agreement_table_list)

# STATS FILTERED POSITIONS ####
## Bejamin Hochberg 
length(Agreement_table_all_chr$Chromosome[Agreement_table_all_chr$Overlap_BH=="TRUE"])/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$ParaMask_fC)/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$NGSparalog_fC_BH)/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$ParaMask_fC)
sum(Agreement_table_all_chr$NGSparalog_fC_BH)

## Bonferroni 
length(Agreement_table_all_chr$Chromosome[Agreement_table_all_chr$Overlap_bonferroni=="TRUE"])/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$ParaMask_fC)/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$NGSparalog_fC_bonferroni)/length(Agreement_table_all_chr$Chromosome)
sum(Agreement_table_all_chr$ParaMask_fC)
sum(Agreement_table_all_chr$NGSparalog_fC_bonferroni)
