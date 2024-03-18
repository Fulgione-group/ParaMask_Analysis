# Create function to compare SVs with ParaMask multicopy-regions
# Additionally, search for SV overlaps with repeats and organella DNA


FNR_ParaMask_function <- function(SV_file,Coordinates_files_list, plastid_list, repeats_list){
  
  # Create a list to store results
  SV_list_name_unprocessed <- list()  
  
  # Loop through chromosomes 1 to 8
  for (chr_processing in chromosome_vector) {
    # Subset data for the current chromosome
    list_in_process <- SV_file[SV_file$Chromosome == chr_processing, ]
    colnames(list_in_process) <- c("Chromosome" ,"ref_start"  ,"ref_end","LENGTH","TYPE")
    
    list_in_process$query_start <- rep(0, length(list_in_process$ref_start))
    list_in_process$query_end <- rep(0, length(list_in_process$ref_end))
    # Save the subset in the list
    Coordinate_file_in_process <- Coordinates_files_list[[chr_processing]]
    
    for (hey in 1:length(Coordinate_file_in_process$query_name)) {
      print(which(list_in_process$ref_start > Coordinate_file_in_process$ref_start[hey] & list_in_process$ref_start < Coordinate_file_in_process$ref_end[hey]))
      list_in_process$query_start[which(list_in_process$ref_start > Coordinate_file_in_process$ref_start[hey] & list_in_process$ref_start < Coordinate_file_in_process$ref_end[hey])] <- Coordinate_file_in_process$query_start[hey]
      list_in_process$query_end[which(list_in_process$ref_start > Coordinate_file_in_process$ref_start[hey] & list_in_process$ref_start < Coordinate_file_in_process$ref_end[hey])] <- Coordinate_file_in_process$query_end[hey]
    }
    colnames(list_in_process) <- c("Chromosome" ,"START"  ,"STOP","LENGTH","TYPE", "query_start", "query_end")
    SV_list_name_unprocessed[[chr_processing]] <- list_in_process
  }
  
  
  
  #Look at all single-copy regions based on the merged set
  SV_list_name_processed <- list()       
  
  for (chr_processing in chromosome_vector) {
    chromosome_in_process <- SV_list_name_unprocessed[[chr_processing]]
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
      chromosome_in_process$Mean_ParaMask_fin_class[pos] <- mean(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]])
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
      chromosome_in_process$number_correct_assigned_SNPs[pos] <-length(which(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==1))
      chromosome_in_process$number_false_assigned_SNPs[pos] <- length(which(EM_ES03_ES04_finalClass[[chr_processing]]$finalClass[EM_ES03_ES04_finalClass[[chr_processing]]$Position > chromosome_in_process$START[pos] & EM_ES03_ES04_finalClass[[chr_processing]]$Position < chromosome_in_process$STOP[pos]]==0))
      
      ### include plastids and repeats ###
      
      plastid_chr_in_process <- plastid_list[[chr_processing]]
      repeat_chr_in_process <- repeats_list[[chr_processing]]
      
      chromosome_in_process$Repeat_class <- rep("NA", nrow(chromosome_in_process))
      chromosome_in_process$Blast_hits <- rep("NA", nrow(chromosome_in_process))
      
      for (plastid_line in 1:length(chromosome_in_process$Chromosome)) {
        
        #Transfer repeat values
        repeat_values <- repeat_chr_in_process$V4[repeat_chr_in_process$V2 >= chromosome_in_process$START[plastid_line] & repeat_chr_in_process$V2 <= chromosome_in_process$STOP[plastid_line]]
        repeat_values <- unique(repeat_values)
        chromosome_in_process$Repeat_class[plastid_line] <- paste(repeat_values, collapse = ", ") 
        
        # Transfer Blast hits
        blast_values <- plastid_chr_in_process$Coverage_classes[plastid_chr_in_process$Window >= chromosome_in_process$query_start[plastid_line] & plastid_chr_in_process$Window <= chromosome_in_process$query_end[plastid_line]]
        chromosome_in_process$Blast_hits[plastid_line] <- paste(blast_values, collapse = ", ")  
      }
      
      SV_list_name_processed[[chr_processing]] <- chromosome_in_process
    }
    
  }
  return(SV_list_name_processed)
}
