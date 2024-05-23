## Prepare input files ####
# READ IN ALL FILES AND STORE THEM IN LIST
# Calculate p-Value and adjust with: 
  # Benjamin Hochburg, p < 0.001
  # Bonferroni, p < 0.05
## ADDITIONALLY FILTER ALL CHROMOSOME FILES FOR OVERLAPPING SNP POSITIONS ONLY 

NGS_paralog_list_original <- list()
EM_ES03_ES04_finalClass_original <- list()


EM_ES03_ES04_finalClass_filtered <- list()
NGS_paralog_list_filtered <- list()

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

# Loop through the chromosomes and assign the corresponding data frames to the list
for (chr in chromosomes) {
  
  #Read all ParaMask results
  PM_currecnt_chr_df <- read.delim(paste("~/PATH/TO/PARAMASK/ES0304_run_finalEMresults.", chr,".finalClass.het", sep = ""))
  EM_ES03_ES04_finalClass_original[[chr]] <- PM_currecnt_chr_df
  
  #Read all NGSParalog results
  NGS_currecnt_chr_df <- read.delim(paste("~/PATH/TO/NGSPARALOG/NGS_out_", chr,".lr", sep=""), header=FALSE)
  colnames(NGS_currecnt_chr_df) <- c("chromosome", "Position", "LR_1", "LR_2", "final_LR")
  # Calculate p-Values and adjust 
  NGS_currecnt_chr_df$pval <- 0.5*pchisq(NGS_currecnt_chr_df$final_LR,df=1,lower.tail=FALSE) 
  NGS_currecnt_chr_df$pval.adj_BH <- p.adjust(NGS_currecnt_chr_df$pval, method="BH") 
  NGS_currecnt_chr_df$pval.adj_bonferroni <- p.adjust(NGS_currecnt_chr_df$pval, method="bonferroni")
  
  # Assing finalClass column
  NGS_currecnt_chr_df$finalClass_BH <- rep(0, length(NGS_currecnt_chr_df$chromosome))
  NGS_currecnt_chr_df$finalClass_bonferroni <- rep(0, length(NGS_currecnt_chr_df$chromosome))
  
  # Assign final class 1 or 0 
  NGS_currecnt_chr_df$finalClass_BH[NGS_currecnt_chr_df$pval.adj_BH < 0.001] <- 1
  NGS_currecnt_chr_df$finalClass_bonferroni[NGS_currecnt_chr_df$pval.adj_bonferroni < 0.05] <- 1
  
  NGS_paralog_list_original[[chr]] <- NGS_currecnt_chr_df
  
  # Filter for SNPs present in ParaMask and ngsparalog 
  overlapping_positions <- intersect(PM_currecnt_chr_df$Position, NGS_currecnt_chr_df$Position)
  
  # Filter the data frames to include only overlapping positions
  filtered_em_df <- PM_currecnt_chr_df[PM_currecnt_chr_df$Position %in% overlapping_positions, ]
  filtered_ngs_df <- NGS_currecnt_chr_df[NGS_currecnt_chr_df$Position %in% overlapping_positions, ]
  
  EM_ES03_ES04_finalClass_filtered[[chr]] <- filtered_em_df
  NGS_paralog_list_filtered[[chr]] <- filtered_ngs_df
}

