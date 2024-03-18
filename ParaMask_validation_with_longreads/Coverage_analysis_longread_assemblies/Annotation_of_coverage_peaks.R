# Assign blast, plastid, SV result to long-read coverage peaks

# Load packages 

# Set working directory
setwd("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/ES04_014/")

# create function to categorize different coverage peaks 
process_coverage_classes <- function(data) {
  cover_classes <- ifelse(data$Mean_Coverage < 20, "coverage_drop",
                          ifelse(data$Mean_Coverage < 100, "mean coverage",
                                 ifelse(data$Mean_Coverage < 200, "2-3-fold coverage", "plastid_DNA")))
  data$Coverage_classes <- cover_classes
  return(data)
}

# Read in generated csv files, including mean coverage in 1500bp windows per chromosome 
tables_ES04 <- lapply(1:8, function(i) {
  data <- read.csv(paste0("ES04_chr", i, "_coverage_mean.csv"))
  process_coverage_classes(data)
})

# Create function which extracts only coverage peaks from coverage data ("2-3-fold coverage", "plastid_DNA")
extract_paralogs <- function(data, chr_num) {
  paralogs <- data[data$Coverage_classes %in% c("2-3-fold coverage", "plastid_DNA"), ]
  colnames(paralogs) <- c("Chr_coverage_self_alignment", "query_start", "Mean_Coverage", "Coverage_classes")
  paralogs$query_end <- NA
  return(paralogs)
}

# Generate list including paralog data for each chromosome-> apply function 
Paralogs_list <- lapply(1:8, function(chr_num) {
  extract_paralogs(tables_ES04[[chr_num]], chr_num)
})

# Generate bed file (size of each potential multicopy region
create_bed_file <- function(paralogs) {
  for (i in 2:nrow(paralogs)) {
    diff_query_start <- paralogs$query_start[i] - paralogs$query_start[i - 1]
    if (diff_query_start < 8000) {
      paralogs$query_end[i] <- NA
    } else {
      new_end_value <- paralogs$query_start[i - 1]
      paralogs$query_end[i] <- new_end_value
    }
  }
  paralogs <- subset(paralogs, !is.na(query_end))
  paralogs$paralog_length <- as.numeric(paralogs$query_end) - as.numeric(paralogs$query_start)
  paralogs <- subset(paralogs, paralog_length != 0)
  paralogs$query_end <- as.numeric(paralogs$query_end)
  paralogs$ref_start <- 0
  paralogs$ref_end <- 0
  return(paralogs)
}

# Apply bed function to each chromosome 
Paralogs_bed_files <- lapply(Paralogs_list, create_bed_file)

# Create function to assign coordinates of reference (match coordinate system file)
assign_reference_coordinates <- function(paralogs, coordinates) {
  for (i in 1:nrow(paralogs)) {
    query_end_value <- paralogs$query_end[i]
    start_index <- which(coordinates$query_start == query_end_value)[1] - 1
    end_index <- which(coordinates$query_start == query_end_value)[1] + 1
    paralogs$ref_start[i] <- coordinates$ref_end[start_index]
    paralogs$ref_end[i] <- coordinates$ref_end[end_index]
  }
  return(paralogs)
}

# Read in files of including coordinate information and apply function 
Coordinates_files_list <- lapply(1:8, function(chr_num) {
  read.delim(paste0("../Files_for_R/Coordinate_files/ES04_014_coordinates_reference_purged.bed"))
})
Paralogs_list_with_ref_coordinates <- Map(assign_reference_coordinates, Paralogs_bed_files, Coordinates_files_list)

#
assign_SVs <- function(paralogs, SV_data) {
  for (i in 1:nrow(paralogs)) {
    query_start <- paralogs$query_start[i]
    query_end <- paralogs$query_end[i]
    overlapping_SVs <- SV_data[SV_data$query_start > query_start & SV_data$query_end < query_end, ]
    paralogs$Coverage_classes[paralogs$Window %in% overlapping_SVs$query_start] <- "verified SV (sniffles/cuteSV)"
  }
  return(paralogs)
}

# Read in SV data and create list for each chromosome 
SV_list <- lapply(1:8, function(chr_num) {
  read.delim(paste0("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/ES04_014/ES04_014_all_merged_SVs_wit_header.txt"))
})

# Read in plastid files and split between chromosomes 
plastid_DNA_files <- lapply(1:8, function(chr_num) {
  read.delim(paste0("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/ES04_014/ES04_014_plastid_DNA.txt"))
})
# Read in rRNA-files and split between chromosomes 
rRNA_files <- lapply(1:8, function(chr_num) {
  read.delim(paste0("~/Studium/Köln_Biological Science/Master Thesis/Bastiaan/ES04_014/ES04_014_rRNA.bed"), header = FALSE)
})


##############
#Assign to coverage files
