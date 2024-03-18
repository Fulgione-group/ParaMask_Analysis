#Load packages
library(GenomeScope)

#create chromosome vector
chromosome_vector <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")


#Read in size of chromosomes 
Chromosome_sizes_ES03_014 <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/VCF_processed_for_R/Single_copy_regions/ES03_014_chromosome_sizes.txt")

Chromosome_sizes_ES03_014_list <- list()

# Loop through chromosomes 1 to 8
for (chr_num in 1:8) {
  # Subset data for the current chromosome
  chr_data <- Chromosome_sizes_ES03_014[Chromosome_sizes_ES03_014$Chromosome == paste0("chr", chr_num), ]
  
  # Save the subset in the list
  Chromosome_sizes_ES03_014_list[[paste0("chr", chr_num)]] <- chr_data
  colnames(Chromosome_sizes_ES03_014_list[[paste0("chr", chr_num)]]) <-  c("Chromosome" ,"ref_start"  ,"ref_end","LENGTH","TYPE")
}



ES03_014_overlapping_SVs_all_types_with_header <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/VCF_processed_for_R/Files_all_SV_types/ES03_014_all_SVs_with_header.txt")
ES03_014_overlapping_SVs_all_types_with_header_filtered <- subset(ES03_014_overlapping_SVs_all_types_with_header, TYPE %in% c("DEL", "DUP", "INV", "INS"))
ES03_014_merg_no_svt_dup_with_header <- ES03_014_overlapping_SVs_all_types_with_header_filtered
ES03_014_single_copy_all_chr_list <- list()
for (chr_num in 1:8) {
  # Subset data for the current chromosome
  chr_data <- ES03_014_merg_no_svt_dup_with_header[ES03_014_merg_no_svt_dup_with_header$Chromosome == paste0("chr", chr_num), ]
  gr <- with(chr_data, GRanges(seqnames = Chromosome, ranges = IRanges(START, STOP)))
  chromosome_length <- max(Chromosome_sizes_ES03_014_list[[paste0("chr", chr_num)]]$ref_end)
  
  full_chromosome <- GRanges(seqnames = paste0("chr", chr_num), ranges = IRanges(1, chromosome_length))
  non_dup_ranges <- setdiff(full_chromosome, gr)
  
  # Extract the coordinates from the non-DUP ranges
  non_dup_coordinates <- as.data.frame(cbind(seqnames(non_dup_ranges), start(non_dup_ranges), end(non_dup_ranges)))
  
  colnames(non_dup_coordinates) <- c("Chromosome" ,"START"  ,"STOP")
  non_dup_coordinates$START <- unlist(non_dup_coordinates$START)
  non_dup_coordinates$STOP <- unlist(non_dup_coordinates$STOP)
  non_dup_coordinates$Chromosome <- paste0("chr", chr_num)
  non_dup_coordinates$LENGTH <- non_dup_coordinates$STOP - non_dup_coordinates$START
  non_dup_coordinates$TYPE <- "single_copy_region"
  
  # Save the subset in the list
  ES03_014_single_copy_all_chr_list[[paste0("chr", chr_num)]] <- non_dup_coordinates
}

ES03_014_single_copy_merged <- bind_rows(ES03_014_single_copy_all_chr_list)


# Apply single copy agreement function using generated single-copy and filtered ParaMask file 
ES03_014_merged_sc_result <- FPR_ParaMask_function(ES03_014_single_copy_merged, ES03_014_EM_finalClass_filtered)
ES03_014_single_copy_merged_processed_data_frame <- dplyr::bind_rows(ES03_014_merged_sc_result)

