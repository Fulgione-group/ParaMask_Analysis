#create function for overlap single-copy regions 

FPR_ParaMask_function <- function(SV_file, EM_ES03_ES04_finalClass){   #Here filter of ES03_014_EM_finalClass_filtered or ES04_014_EM_finalClass_filtered

  SV_list_name_processed <- list()
  
  
  for (chr_processing in chromosome_vector) {
    chromosome_in_process <- SV_file[SV_file$Chromosome == chr_processing, ]
    chromosome_in_process$Mean_ParaMask_fin_class <- rep(0, length(chromosome_in_process$Chromosome))
    #chromosome_in_process$Mean_ParaMask_EM_class <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$Mean_ParaMask_coverage <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$Mean_ParaMask_coverage_het_seeds <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$Mean_ParaMask_coverage_het_uncertain<- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$sum_SNP_seed <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$sum_SNP_uncertain <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$sum_SNP_single_copy <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$proportion_SNP_seed <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$proportion_SNP_uncertain <- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$proportion_SNP_single_copy <- rep(0, length(chromosome_in_process$Chromosome))
    
    chromosome_in_process$number_all_SNPs<- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$number_correct_assigned_SNPs<- rep(0, length(chromosome_in_process$Chromosome))
    chromosome_in_process$number_false_assigned_SNPs <- rep(0, length(chromosome_in_process$Chromosome))
    
    ####### compare
    
    
    for (pos in 1:length(chromosome_in_process$Chromosome)) {
      print(mean(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]))
      
      chromosome_in_process$Mean_ParaMask_fin_class[pos] <- 1-mean(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]])
      #chromosome_in_process$Mean_ParaMask_EM_class[pos] <- mean(EM_ES03_ES04_finalClass[[chr_processing]]$EM_class[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]])
      
      chromosome_in_process$Mean_ParaMask_coverage[pos] <- mean(EM_ES03_ES04_finalClass[[chr_processing]]$Mean.coverage[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]])
      
      chromosome_in_process$Mean_ParaMask_coverage_het_seeds[pos] <- mean(EM_ES03_ES04_finalClass[[chr_processing]][
        EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] &
          EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos],]$Mean.coverage.het[EM_ES03_ES04_finalClass[[chr_processing]][
            EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] &
              EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos],]$EM_class==2])
      
      chromosome_in_process$Mean_ParaMask_coverage_het_uncertain[pos] <- mean(EM_ES03_ES04_finalClass[[chr_processing]][
        EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] &
          EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos],]$Mean.coverage.het[EM_ES03_ES04_finalClass[[chr_processing]][
            EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] &
              EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos],]$EM_class==1])
      
      
      chromosome_in_process$sum_SNP_seed[pos] <-  length(which(EM_ES03_ES04_finalClass[[chr_processing]]$EM_class[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==2))
      chromosome_in_process$sum_SNP_uncertain[pos] <-length(which(EM_ES03_ES04_finalClass[[chr_processing]]$EM_class[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==1))
      chromosome_in_process$sum_SNP_single_copy[pos] <-length(which(EM_ES03_ES04_finalClass[[chr_processing]]$EM_class[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==0))
      
      chromosome_in_process$proportion_SNP_seed[pos] <-   chromosome_in_process$sum_SNP_seed[pos] / (chromosome_in_process$sum_SNP_seed[pos]+ chromosome_in_process$sum_SNP_uncertain[pos]+chromosome_in_process$sum_SNP_single_copy[pos])
      chromosome_in_process$proportion_SNP_uncertain[pos] <-   chromosome_in_process$sum_SNP_uncertain[pos] / (chromosome_in_process$sum_SNP_seed[pos]+ chromosome_in_process$sum_SNP_uncertain[pos]+chromosome_in_process$sum_SNP_single_copy[pos])
      chromosome_in_process$proportion_SNP_single_copy[pos] <-   chromosome_in_process$sum_SNP_single_copy[pos] / (chromosome_in_process$sum_SNP_seed[pos]+ chromosome_in_process$sum_SNP_uncertain[pos]+chromosome_in_process$sum_SNP_single_copy[pos])
      
      chromosome_in_process$number_all_SNPs[pos] <- length(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]])
      chromosome_in_process$number_correct_assigned_SNPs[pos] <-length(which(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==0))
      chromosome_in_process$number_false_assigned_SNPs[pos] <- length(which(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==1))
      
      SV_list_name_processed[[chr_processing]] <- chromosome_in_process
    }
    
  }
  return(SV_list_name_processed)
}
