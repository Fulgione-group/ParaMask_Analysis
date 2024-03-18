#Run analysis of agreement in multicopy regions

#defined function FNR_function

#Read in Files: SV_file, coordinate file, repeat_file, organella_file, ParaMask het file
ES03_014_SV_loop_merged <- read.delim("VCF_processed_for_R/Files_only_with_duplications/ES03_014_only_overlapping_dup_with_header.txt")

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

#Read HET files
EM_ES03_ES04_chr1.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr1.finalClass.het")
EM_ES03_ES04_chr2.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr2.finalClass.het")
EM_ES03_ES04_chr3.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr3.finalClass.het")
EM_ES03_ES04_chr4.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr4.finalClass.het")
EM_ES03_ES04_chr5.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr5.finalClass.het")
EM_ES03_ES04_chr6.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr6.finalClass.het")
EM_ES03_ES04_chr7.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr7.finalClass.het")
EM_ES03_ES04_chr8.finalClass <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/ES0304_run_finalEMresults.chr8.finalClass.het")


EM_ES03_ES04_finalClass <- list()
# Loop through the chromosomes and assign the corresponding data frames to the list
for (chr in chromosomes) {
  df_name <- paste("EM_ES03_ES04_", chr, ".finalClass", sep = "")
  EM_ES03_ES04_finalClass[[chr]] <- get(df_name)
}
#### Prepare Data with coordinates and plastids/repeats ####

### Read in coordinate files to match self-alignment/repeat annotation coordinates with reference coordinate #####
chromosome_vector <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")
setwd("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/")

### for ES03-014
ES03_014_Coordinates_reference_scaffolds_with_headers <- read.delim("Coordinate_files/ES03_014_coordinates_reference_unpurged.bed")

ES03_014_Coordinates_files_list <- list()


# Loop through chromosomes 1 to 8
for (chr_num in 1:8) {
  # Subset data for the current chromosome
  chr_data <- ES03_014_Coordinates_reference_scaffolds_with_headers[ES03_014_Coordinates_reference_scaffolds_with_headers$query_name == paste0("chr", chr_num), ]
  colnames(chr_data) <- c("query_name","ref_start","ref_end","reference","query_start","query_end")
  chr_data <- chr_data[order(chr_data$query_start), ]
  # Save the subset in the list
  ES03_014_Coordinates_files_list[[paste0("chr", chr_num)]] <- chr_data
}


#####include all repeats and plastids#######
ES03_014_repeats <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Repeat_files/Repeats_ES03_with_reference_coordinates.bed", quote="\"", comment.char="")
ES03_014_repeats <- na.omit(ES03_014_repeats)
ES03_014_plastid_DNA <- read.delim("Blast_search_files/ES03_014_plastid_DNA.txt")
ES03_014_plastid_list <- list()
ES03_014_repeats_list <- list()

for (chr_processing in chromosome_vector) {
  repeat_in_process <- ES03_014_repeats[ES03_014_repeats$V1==chr_processing,]
  ES03_014_repeats_list[[chr_processing]] <- repeat_in_process
  plastid_chromosome <- paste(chr_processing,"_RagTag", sep = "")
  plastid_in_process <- ES03_014_plastid_DNA[ES03_014_plastid_DNA$Chr== paste(chr_processing,"_RagTag", sep = ""),]
  ES03_014_plastid_list[[chr_processing]] <- plastid_in_process
}

### for ES04-014
ES04_014_repeats <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Repeat_files/Repeats_ES04_with_reference_coordinates.bed", quote="\"", comment.char="")
ES04_014_repeats <- na.omit(ES04_014_repeats)
ES04_014_plastid_DNA <- read.delim("Blast_search_files/ES04_014_plastid_DNA.txt")

ES04_014_plastid_list <- list()
ES04_014_repeats_list <- list()

for (chr_processing in chromosome_vector) {
  repeat_in_process <- ES04_014_repeats[ES04_014_repeats$V1==chr_processing,]
  ES04_014_repeats_list[[chr_processing]] <- repeat_in_process
  plastid_chromosome <- paste(chr_processing,"_RagTag", sep = "")
  plastid_in_process <- ES04_014_plastid_DNA[ES04_014_plastid_DNA$Chr==plastid_chromosome,]
  ES04_014_plastid_list[[chr_processing]] <- plastid_in_process
}



##### Read in all SV sets of different tools #####
ES04_014_SV_loop_merged <- read.delim("VCF_processed_for_R/Files_only_with_duplications/ES04_014_only_overlapping_dup_with_header.txt")




#Run function using all files provided 
ES03_014_merged_result <- FNR_ParaMask_function(ES03_014_SV_loop_merged, ES03_014_Coordinates_files_list,ES03_014_plastid_list, ES03_014_repeats_list)
ES03_014_SV_loop_merged_processed_data_frame <- dplyr::bind_rows(ES03_014_merged_result)

ES04_014_merged_result <- FNR_ParaMask_function(ES04_014_SV_loop_merged, ES04_014_Coordinates_files_list,ES04_014_plastid_list, ES04_014_repeats_list)
ES04_014_SV_loop_merged_processed_data_frame <- dplyr::bind_rows(ES04_014_merged_result)


#Process data frame further
ES03_014_SV_loop_merged_processed_data_frame$SV_tool <- rep("merged", length(ES03_014_SV_loop_merged_processed_data_frame$Chromosome))
ES04_014_SV_loop_merged_processed_data_frame$SV_tool <- rep("merged", length(ES04_014_SV_loop_merged_processed_data_frame$Chromosome))

ES03_014_SV_loop_merged_processed_data_frame$accession <- rep("ES03-014", length(ES03_014_SV_loop_merged_processed_data_frame$Chromosome))
ES04_014_SV_loop_merged_processed_data_frame$accession <- rep("ES04-014", length(ES04_014_SV_loop_merged_processed_data_frame$Chromosome))


#Bind tables
Both_accessions_merged_SVs_df <- rbind(ES03_014_SV_loop_merged_processed_data_frame, ES04_014_SV_loop_merged_processed_data_frame)
write.table(Both_accessions_merged_SVs_df, "~/Studium/Köln_Biological Science/Master Thesis/ParaMask/Results/Paralog_agreement_validation_merged_SVs.txt")
