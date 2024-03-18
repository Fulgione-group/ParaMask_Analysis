#What is ParaMask calling? -> Overlap between multicopy SNPs and Repeats/SVs
#load packages 
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)

cbPalette_9 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")

##### 1) Filter the HET for ES03 seeds only and minimum single-copy length of 35bp #####
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")


# Define the base file path
base_path <- "~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Final/"

# Initialize lists to store dataframes
EM_ES03_ES04_finalClass_original <- list()
Clusters_final_ES03_as_seed_list <- list()

# Loop through chromosomes and read files
for (chr in chromosomes) {
  # Read HET files
  het_file <- paste0(base_path, "ES0304_run_finalEMresults.", chr, ".finalClass.het")
  EM_ES03_ES04_finalClass_original[[chr]] <- read.delim(het_file)
  
  # Read cluster files
  cluster_file <- paste0(base_path, "cluster_files_filtered/Clusters_ES03_as_seed/clusters_final_ES03_as_seed_", chr, ".txt")
  Clusters_final_ES03_as_seed_list[[chr]] <- read.delim(cluster_file, header = FALSE)
}


ES03_014_EM_finalClass_only_ES03_seed <- EM_ES03_ES04_finalClass_original
# Loop through the chromosomes and assign the corresponding data frames to the list
for (chr in chromosomes) {
  subset_het_filtered <- EM_ES03_ES04_finalClass_original[[chr]]
  subset_cluster_filtered <- Clusters_final_ES03_as_seed_list[[chr]]
  unique(subset_cluster_filtered$V2)
  
  #read inbed file for size filter [35bp]
  file_path_bed_files <- paste0(base_path_bed_files, chr, ".finalClass.bed")
  current_bed_file <- read.delim(file_path_bed_files, header = TRUE)
  
  # Extract the rows in which are paralogs present
  current_paralog_df <- current_bed_file[current_bed_file$type.0.single.copy.1.multi.copy == 1, ]
  
  #calculate length
  current_paralog_df$length <- current_paralog_df$End-current_paralog_df$Start
  size_filter_df <- current_paralog_df[which(current_paralog_df$length<35),]    #filter for size
  length(size_filter_df$Chromosome)
  length(subset_cluster_filtered$V1)
  subset_cluster_filtered_size <- subset_cluster_filtered[-which(unique(subset_cluster_filtered$V3) %in% size_filter_df$Cluster), ]
  
  #Set all final class values to 0 
  subset_het_filtered$finalClass <- 0
  
  #Only assign ES03 seed cluster as paralogous SNP
  subset_het_filtered$finalClass[which(subset_het_filtered$Position%in%subset_cluster_filtered_size$V2)] <- 1
  subset_het_filtered$SV_class <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Repeat_class <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Blast_hit <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Repeat_priority <- rep(NA, length(subset_het_filtered$Chromosome))
  ES03_014_EM_finalClass_only_ES03_seed[[chr]] <- subset_het_filtered
  Clusters_final_ES03_as_seed_list[[chr]] <- subset_cluster_filtered_size
}

ES03_014_EM_finalClass_only_ES03_seed_data_frame <- bind_rows(ES03_014_EM_finalClass_only_ES03_seed)
EM_ES03_ES04_finalClass_original_data_frame <- bind_rows(EM_ES03_ES04_finalClass_original)

#Count length of paralogous SNPs before and after filtering 
length(ES03_014_EM_finalClass_only_ES03_seed_data_frame$Chromosome[ES03_014_EM_finalClass_only_ES03_seed_data_frame$finalClass==1])
length(EM_ES03_ES04_finalClass_original_data_frame$Chromosome[EM_ES03_ES04_finalClass_original_data_frame$finalClass==1])




### Filtered HET File: ES03_014_EM_finalClass_filtered; 
##### 2) Read in all necessary files and split for each chr #####
#### Coordinates for assembly files 
#### SV calls 
#### Repeats
#### Blast hits

# Read in SV files
ES03_014_overlapping_SVs_all_types_with_header <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/VCF_processed_for_R/Files_all_SV_types/ES03_014_only_overlapping_SVs_with_header.txt")
ES03_014_overlapping_SVs_all_types_with_header <- ES03_014_overlapping_SVs_all_types_with_header[ES03_014_overlapping_SVs_all_types_with_header$TYPE!="TRA",]

ES03_014_overlapping_SVs_all_types_with_header$LENGTH <- abs(as.integer(ES03_014_overlapping_SVs_all_types_with_header$LENGTH))

#Read in Coordinate Files
ES03_014_Coordinates_reference_scaffolds_with_headers <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Coordinate_files/ES03_014_coordinates_reference_unpurged.bed")

#Read in Repeats
Arabis_repeats <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/Files_for_R/Repeat_files/Repeats_ES04_with_reference_coordinates.bed", quote="\"", comment.char="")
Arabis_repeats <- na.omit(Arabis_repeats)

#Read in Plastid Files (not used in the anaylsis)
Arabis_plastid_DNA <- read.delim("Blast_search_files/ES03_014_plastid_DNA.txt")

#Create list for all chromosomes of each category
ES03_014_overlapping_SVs_all_types_list <- list()
ES03_014_Coordinates_files_list <- list()
Arabis_plastid_list <- list()
Arabis_repeats_list <- list() 

for (chr_processing in chromosomes) {
  #SV calls
  SV_data <- ES03_014_overlapping_SVs_all_types_with_header[ES03_014_overlapping_SVs_all_types_with_header$Chromosome== chr_processing, ]
  # Save the subset in the list
  SV_data$ParaMask_overlap <- rep("No", length(SV_data$Chromosome))
  ES03_014_overlapping_SVs_all_types_list[[chr_processing]] <- SV_data
  
  
  #Coordinates
  chr_data <- ES03_014_Coordinates_reference_scaffolds_with_headers[ES03_014_Coordinates_reference_scaffolds_with_headers$query_name == chr_processing, ]
  colnames(chr_data) <- c("query_name","ref_start","ref_end","reference","query_start","query_end")
  chr_data <- chr_data[order(chr_data$query_start), ]
  ES03_014_Coordinates_files_list[[chr_processing]] <- chr_data
  
  #Repeats
  repeat_in_process <- Arabis_repeats[Arabis_repeats$V1==chr_processing,]
  Arabis_repeats_list[[chr_processing]] <- repeat_in_process
}

##### 3) Assign SV, Repeat, Blast hits for each paralogous SNP (HET file filtered as input): #####

### 3.1) For each paralogous SNP ###
#if you need filtered HET file for other analysis since 'ES03_014_EM_finalClass_only_ES03_seed' will be updated for only paralogous SNPs
ES03_014_EM_finalClass_filtered <- ES03_014_EM_finalClass_only_ES03_seed

#start adding all necessary columns and remove all not paralogous SNPs from ES03_014_EM_finalClass_only_ES03_seed
for (chr in chromosomes) {
  subset_het <- ES03_014_EM_finalClass_filtered[[chr]]
  subset_het_filtered <- subset_het[subset_het$finalClass==1,]
  subset_het_filtered$SV_class <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Repeat_class <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Blast_hit <- rep(NA, length(subset_het_filtered$Chromosome))
  subset_het_filtered$Repeat_priority <- rep(NA, length(subset_het_filtered$Chromosome))
  
  ES03_014_EM_finalClass_only_ES03_seed[[chr]] <- subset_het_filtered
}

###Now start analysing the overlap of Repeats/SVs with paralogous SNPs###

HET_files_list_processed_seed_filtered <- list()
for (chr in chromosomes) {
 print(chr)
  subset_paramask <- ES03_014_EM_finalClass_only_ES03_seed[[chr]]
  subset_repeats <- Arabis_repeats_list[[chr]]
  subset_blast <- Arabis_plastid_list[[chr]]
  for (position in 1:length(ES03_014_overlapping_SVs_all_types_list[[chr]]$TYPE)) {
    start_condition <-  ES03_014_overlapping_SVs_all_types_list[[chr]]$START[position] 
    stop_condition <-  ES03_014_overlapping_SVs_all_types_list[[chr]]$STOP[position]
    subset_paramask$SV_class[subset_paramask$Position > start_condition & subset_paramask$Position  < stop_condition] <- ES03_014_overlapping_SVs_all_types_list[[chr]]$TYPE[position]
    subset_paramask$Repeat_priority[subset_paramask$Position > start_condition & subset_paramask$Position  < stop_condition] <- ES03_014_overlapping_SVs_all_types_list[[chr]]$TYPE[position]
    
  }
  for (position in 1:length(subset_repeats$V1)) {
    start_condition_repeat <-  subset_repeats$V2[position] 
    stop_condition_repeat <-  subset_repeats$V3[position]
    subset_paramask$Repeat_class[subset_paramask$Position > start_condition_repeat & subset_paramask$Position  < stop_condition_repeat] <- subset_repeats$V4[position]
    subset_paramask$Repeat_priority[subset_paramask$Position > start_condition_repeat & subset_paramask$Position  < stop_condition_repeat] <- subset_repeats$V4[position]
  }
  HET_files_list_processed_seed_filtered[[chr]] <- subset_paramask
}
ParaMask_SNP_comparison <- bind_rows(HET_files_list_processed_seed_filtered)




selected_colnames <- c("Chromosome", "Position", "EM_class", "finalClass", "SV_class", "Repeat_class", "Repeat_priority")
ParaMask_SNP_subset <- ParaMask_SNP_comparison[selected_colnames]

table(ParaMask_SNP_comparison$SV_class)
table(ParaMask_SNP_comparison$Repeat_class)
table(ParaMask_SNP_comparison$Repeat_priority)
unique(ParaMask_SNP_comparison$Repeat_priority)

# Specify the desired column names


# Your provided list of transposable element categories
transposable_elements <- c("DNA/MULE-MuDR", "DNA/TcMar-Stowaway", "DNA/hAT-Ac", "LTR/Gypsy", "DUP", "DEL", "LINE/L1", "RC/Helitron", 
                           "DNA/CMC-EnSpm", "LTR/Copia", "DNA/TcMar-Pogo", "TRA", "DNA", "DNA/PIF-Harbinger", "LTR/Cassandra", "DNA/hAT-Tag1", 
                           "DNA/hAT-Tip100", "SINE/tRNA", "LTR/Pao", "SINE", "INV", "DNA/hAT-Charlie", "LINE/L1-Tx1", "LTR", "Retroposon/L1-dep", 
                           "Retroposon", "LTR/Caulimovirus", "DNA/hAT", "DNA/Maverick", "DNA/TcMar-Fot1", "rRNA", "LTR/ERVK", 
                           "LINE", "tRNA", "LINE/R1", "Satellite")

# Categorize into classes
# Categorize into classes
retrotransposons <- c("LTR/Copia", "LTR/Gypsy", "LTR/ERVK", "LINE/L1", "LTR/Pao", "LINE/L1-Tx1", "LINE/R1", "LTR", "LTR/Caulimovirus","LTR/Cassandra", "Retroposon/L1-dep" ,"Retroposon" , "LTR/Caulimovirus", "SINE/tRNA","SINE","LINE" )
dna_transposons <- c("RC/Helitron","DNA/MULE-MuDR", "DNA/TcMar-Stowaway", "DNA/hAT-Ac", "DNA/CMC-EnSpm", "DNA/TcMar-Pogo", "DNA/PIF-Harbinger", "DNA/hAT-Tag1", "DNA/hAT-Tip100", "DNA/hAT-Charlie", "DNA/Maverick", "DNA/TcMar-Fot1", "DNA/hAT")
deletion <- "DEL"
duplication <- "DUP"
inversion <- "INV"
translocation <- "TRA"
other_elements <- setdiff(transposable_elements, c(retrotransposons, dna_transposons, deletion, duplication, inversion, translocation))


# Display the categorization
cat("Retrotransposons:", retrotransposons, "\n\n")
cat("DNA transposons:", dna_transposons, "\n\n")
cat("Other elements:", other_elements, "\n")



# Assuming your data frame is named ParaMask_SNP_comparison
# and Repeat_priority is the column containing transposable element classes

ParaMask_SNP_subset$Transposable_Class <- "NA"

# Assign categories to the new column
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% retrotransposons] <- "Retrotransposon"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% dna_transposons] <- "DNA Transposon"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% Deletion] <- "Deletion"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% Duplication] <- "Duplication"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% Inversion] <- "Inversion"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% other_elements] <- "Other Element"
ParaMask_SNP_subset$Transposable_Class[ParaMask_SNP_subset$Repeat_priority %in% Translocation] <- "Translocation"
unique(ParaMask_SNP_subset$Transposable_Class)


# Save dataframe with translocations included
write.table(ParaMask_SNP_subset, "what_is_paramask_calling_original.txt")

sum(length(ParaMask_SNP_subset$Chromosome[ParaMask_SNP_subset$Transposable_Class=="NA"]))/length(ParaMask_SNP_subset$Transposable_Class)

# Remove rows with "Translocation" and "NA" from the dataframe
ParaMask_SNP_subset_ohne_translocation <- ParaMask_SNP_subset[ParaMask_SNP_subset$Transposable_Class != "Translocation", ]
ParaMask_SNP_subset_ohne_trans_NA <- ParaMask_SNP_subset_ohne_translocation[ParaMask_SNP_subset_ohne_translocation$Transposable_Class != "NA", ]



write.table(ParaMask_SNP_subset_ohne_trans_NA, "what_is_paramask_calling_wo_TRA_and_NA.txt")
