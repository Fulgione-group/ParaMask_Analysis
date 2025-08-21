
##wrapper to write fastas
wrap_sequence <- function(seq_chars, width = 60) {
  seq_string <- paste(seq_chars, collapse = "")
  n <- nchar(seq_string)
  starts <- seq(1, n, by = width)
  ends <- pmin(starts + width - 1, n)
  return(substring(seq_string, starts, ends))
}
# Set output directory
output_dir <- "$PATH_TO_OUT_DIR/fasta_output"

bases<-c("A", "T","C", "G")
probs <- c(0.3,0.3, 0.2,0.2)
# create a reference sequence
seq<- sample(x = bases, replace = T, prob = probs, size = 1000000)
header <- ">reference"
outfile <- file.path(output_dir,"reference.fasta")
wrapped_seq <- wrap_sequence(seq)
writeLines(c(header, wrapped_seq), con = outfile)





#create a derived seqeunce
seq_derived<- sample(x = bases, replace = T, prob = probs, size = 1000000)
#derived seqeunce must be different
while(sum(seq==seq_derived)>0){
  seq_derived[seq==seq_derived] <- sample(x = bases, replace = T, prob = probs, size = sum(seq==seq_derived))
}
sum(seq==seq_derived)
header <- ">derived"
outfile <- file.path(output_dir,"derived.fasta")
wrapped_seq <- wrap_sequence(seq)
writeLines(c(header, wrapped_seq), con = outfile)


####for HW

###
sum(l_numbers_SC)


#create sequences, copy the sequence according to position and length of each chunk and also look for mutations and copy from the derived sequence


#for a reduced sampling set 
subsamplesize <-15

for(k in 1:3){

  #read all duplication files 0 and 2
  files_SV_0<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/", pattern = paste0(".*_0_.*SV.*rep", k,".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SV_0))
  files_SV_0_sorted <- files_SV_0[order(b_numbers)]
  l_numbers_SV <- as.numeric(sub(".*_l(\\d+)_.*", "\\1", files_SV_0_sorted))
  sum(l_numbers_SV)

  files_SV_2<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/", pattern = paste0(".*_2_.*SV.*rep", k,".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SV_2))
  files_SV_2_sorted <- files_SV_2[order(b_numbers)]

  #read all single-copy files
  files_SC<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/", pattern = paste0("*SC.*rep", k,".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SC))
  files_SC_sorted <- files_SC[order(b_numbers)]
  l_numbers_SC <- as.numeric(sub(".*_l(\\d+)_.*", "\\1", files_SC_sorted))
  sum(l_numbers_SC)

  positions_snps <- c()
  positions_snps_wd <- c()
  fastas <-list()
  i<-1
  SC <- read.table(file = files_SC_sorted[i], header = F, sep = "")
  n_haps<- ncol(SC)
  subsample<- sample(x = ((n_haps-1)/2), replace = F, size = subsamplesize)
  subsample <- subsample*2
  subsample <- sort(c(subsample, subsample-1))+1
  for(j in 2:n_haps){
    fastas[[(j-1)]]  <- c(paste0("> Rep",as.character(k) ,"_Sample", ceiling((j-1)/2),"_Hap", (2-((j-1)%%2))))
  }

  pos_counter <- 1
  pos_counter_wd <- 1
  i<-1
  for(i in 1:length(files_SV_0)){
    SC <- read.table(file = files_SC_sorted[i], header = F, sep = "")
    if(ncol(SC)>1 & nrow(SC)>1){
      positions_snps <- rbind(positions_snps, cbind(SC$V1+pos_counter, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]), "single-copy"))
      positions_snps_wd <- rbind(positions_snps_wd, cbind(SC$V1+pos_counter_wd, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]),"single-copy"))
    }
    l_SC <- l_numbers_SC[i]
    if(l_SC>0){
      for(j in 2:n_haps){
        new_seq <- seq[pos_counter:(pos_counter+l_SC-1)]
        if(ncol(SC)>1){
          derived_pos <- SC$V1[SC[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      pos_counter <- pos_counter + l_SC
    }
    pos_counter_wd <- pos_counter_wd + l_SC
    SV_0 <- read.table(file = files_SV_0_sorted[i], header = F, sep = "")
    l_SV <- l_numbers_SV[i]
    if(l_SV>0){
      if(ncol(SV_0)>1 & nrow(SV_0)>1){
        positions_snps <- rbind(positions_snps, cbind(SV_0$V1+pos_counter, rowSums(SV_0[, 2:n_haps]), rowSums(SV_0[, subsample]), "multi-copy"))
        positions_snps_wd <- rbind(positions_snps_wd, cbind(SV_0$V1+pos_counter_wd, rowSums(SV_0[, 2:n_haps]),rowSums(SV_0[, subsample]), "multi-copy"))
      }
      for(j in 2:n_haps){
        new_seq<- seq[pos_counter:(pos_counter+l_SV-1)]
        if(ncol(SV_0)>1){
          derived_pos<-SV_0$V1[SV_0[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      pos_counter_wd <- pos_counter_wd + l_SV
      SV_2 <- read.table(file = files_SV_2_sorted[i], header = F, sep = "")
      if(ncol(SV_2)>1 & nrow(SV_2)>1){
        positions_snps <- rbind(positions_snps, cbind(SV_2$V1+pos_counter, rowSums(SV_2[, 2:n_haps]), rowSums(SV_2[, subsample]), "multi-copy"))
        positions_snps_wd <- rbind(positions_snps_wd, cbind(SV_2$V1+pos_counter_wd, rowSums(SV_2[, 2:n_haps]),rowSums(SV_2[, subsample]), "multi-copy"))
      }
      for(j in 2:n_haps){
        new_seq<- seq[pos_counter:(pos_counter+l_SV-1)]
        if(ncol(SV_2)>1){
          derived_pos<-SV_2$V1[SV_2[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      pos_counter <- pos_counter +l_SV
      pos_counter_wd <- pos_counter_wd + l_SV
    }
    print(i)
  }

  SC <- read.table(file = files_SC_sorted[length(files_SC_sorted)], header = F, sep = "")
  if(ncol(SC)>1 & nrow(SC)>1){
    positions_snps <- rbind(positions_snps, cbind(SC$V1+pos_counter, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]), "single-copy"))
    positions_snps_wd <- rbind(positions_snps_wd, cbind(SC$V1+pos_counter_wd, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]),"single-copy"))
  }
  l_SC <- l_numbers_SC[length(files_SC_sorted)]
  for(j in 2:n_haps){
    new_seq <- seq[pos_counter:(pos_counter+l_SC-1)]
    if(ncol(SC)>1){
      derived_pos <- SC$V1[SC[j]==1]+1
      new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
    }
    fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
  }

  for(i in 1:length(fastas)){
    fastas[[i]] <- fastas[[i]][!is.na(fastas[[i]])]
  }
  positions_snps <- as.data.frame(positions_snps)
  positions_snps$V1 <- as.numeric(positions_snps$V1)
  positions_snps <- positions_snps[order(positions_snps$V1),]
  # positions_snps <- unique(positions_snps)

  write.table(positions_snps, file = paste0(output_dir, "/Position_files_HW","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)

  positions_snps_wd <- as.data.frame(positions_snps_wd)
  positions_snps_wd$V1 <- as.numeric(positions_snps_wd$V1)
  positions_snps_wd <- positions_snps_wd[order(positions_snps_wd$V1),]
  # positions_snps_wd <- unique(positions_snps_wd)

  write.table(positions_snps_wd, file = paste0(output_dir, "/Position_wd_files_HW","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(subsample-1, file = paste0(output_dir, "/Subsampled_haplotypes_HW","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)

  # Function to wrap sequence into lines of 60 characters


  # Loop over all fastas
  i <-1
  for (i in seq_along(fastas)) {
    header <- fastas[[i]][1]
    sequence <- fastas[[i]][-1]

    # Wrap the sequence
    wrapped_seq <- wrap_sequence(sequence)

    # Define output filename (e.g., Sample1.fasta)
    outfile <- file.path(output_dir, paste0("HW","_Rep",k,"_Sample", ceiling(i/2), "_hap", (2-((i)%%2)) , ".fasta"))

    # Write to file
    writeLines(c(header, wrapped_seq), con = outfile)
  }
}



# for Fis=0.9



###
sum(l_numbers_SC)




Fis <- 0.9
#for a reduced sampling set
subsamplesize <-15
for(k in 1:3){
  #read all duplication files 0 and 2
  files_SV_0<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/Fis0.9/", pattern = paste0(".*_0_.*SV.*rep",k , ".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SV_0))
  files_SV_0_sorted <- files_SV_0[order(b_numbers)]
  l_numbers_SV <- as.numeric(sub(".*_l(\\d+)_.*", "\\1", files_SV_0_sorted))
  sum(l_numbers_SV)

  files_SV_2<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/Fis0.9/", pattern = paste0(".*_2_.*SV.*rep",k , ".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SV_2))
  files_SV_2_sorted <- files_SV_2[order(b_numbers)]



  #read all single-copy files
  files_SC<-list.files(path = "$PATH_TO_SEDUS_OUTPUT/Fis0.9/", pattern = paste0(".*SC.*rep",k,".*transformed\\.txt"), full.names = T)
  # Extract the b<number> part and convert to numeric
  b_numbers <- as.numeric(sub(".*_b(\\d+)_.*", "\\1", files_SC))
  files_SC_sorted <- files_SC[order(b_numbers)]
  l_numbers_SC <- as.numeric(sub(".*_l(\\d+)_.*", "\\1", files_SC_sorted))
  sum(l_numbers_SC)

  positions_snps <- c()
  positions_snps_wd <- c()
  fastas <-list()
  i<-1
  SC <- read.table(file = files_SC_sorted[i], header = F, sep = "")
  n_haps<- ncol(SC)
  subsample<- sample(x = ((n_haps-1)/2), replace = F, size = subsamplesize)
  subsample <- subsample*2
  subsample <- sort(c(subsample, subsample-1))+1
  for(j in 2:n_haps){
    fastas[[(j-1)]]  <- c(paste0("> Rep",as.character(k) ,"_Sample", ceiling((j-1)/2),"_Hap", (2-((j-1)%%2))))
  }
  pos_counter <- 1
  pos_counter_wd <- 1
  i<-1
  for(i in 1:length(files_SV_0)){
    SC <- read.table(file = files_SC_sorted[i], header = F, sep = "")
    l_SC <- l_numbers_SC[i]
    if(l_SC>0){
      if(ncol(SC)>1){
        for(m in 1:nrow(SC)){
          for(n in 1:((n_haps-1)/2)){
            if(sum(SC[m, (n*2):(n*2+1)])==1){
              hetorhom <- sample(x = c(1,2), size = 1, replace = T, prob = c((1-Fis), Fis))
              if(hetorhom==2){
                whichhom <- sample(x = c(0,1), size = 1, replace = T, prob = c(0.5,0.5))
                SC[m, (n*2):(n*2+1)] <- c(whichhom,whichhom)
              }
            }
          }
        }
        SC <- SC[rowSums(SC[, 2:n_haps]) > 0,]
        if(ncol(SC)>1 & nrow(SC)>1){
          positions_snps <- rbind(positions_snps, cbind(SC$V1+pos_counter, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]),"single-copy"))
          positions_snps_wd <- rbind(positions_snps_wd, cbind(SC$V1+pos_counter_wd, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]),"single-copy"))
        }
      }
      for(j in 2:n_haps){
        new_seq <- seq[pos_counter:(pos_counter+l_SC-1)]
        if(ncol(SC)>1){
          derived_pos <- SC$V1[SC[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      pos_counter <- pos_counter + l_SC
      pos_counter_wd <- pos_counter_wd + l_SC
    }
    SV_0 <- read.table(file = files_SV_0_sorted[i], header = F, sep = "")
    l_SV <- l_numbers_SV[i]
    if(l_SV>0){
      if(ncol(SV_0)>1){
        for(m in 1:nrow(SV_0)){
          for(n in 1:((n_haps-1)/2)){
            if(sum(SV_0[m, (n*2):(n*2+1)])==1){
              hetorhom <- sample(x = c(1,2), size = 1, replace = T, prob = c((1-Fis), Fis))
              if(hetorhom==2){
                whichhom <- sample(x = c(0,1), size = 1, replace = T, prob = c(0.5,0.5))
                SV_0[m, (n*2):(n*2+1)] <- c(whichhom,whichhom)
              }
            }
          }
        }
        SV_0 <- SV_0[rowSums(SV_0[, 2:n_haps]) > 0, ]
        if(ncol(SV_0)>1 & nrow(SV_0)>1){
          positions_snps <- rbind(positions_snps, cbind(SV_0$V1+pos_counter, rowSums(SV_0[, 2:n_haps]) , rowSums(SV_0[, subsample]), "multi-copy"))
          positions_snps_wd <- rbind(positions_snps_wd, cbind(SV_0$V1+pos_counter_wd, rowSums(SV_0[, 2:n_haps]), rowSums(SV_0[, subsample]), "multi-copy"))
        }
      }
      pos_counter_wd <- pos_counter_wd + l_SV
      for(j in 2:n_haps){
        new_seq<- seq[pos_counter:(pos_counter+l_SV-1)]
        if(ncol(SV_0)>1){
          derived_pos<-SV_0$V1[SV_0[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      SV_2 <- read.table(file = files_SV_2_sorted[i], header = F, sep = "")
      if(ncol(SV_2)>1){
        for(m in 1:nrow(SV_2)){
          for(n in 1:((n_haps-1)/2)){
            if(sum(SV_2[m, (n*2):(n*2+1)])==1){
              hetorhom <- sample(x = c(1,2), size = 1, replace = T, prob = c((1-Fis), Fis))
              if(hetorhom==2){
                whichhom <- sample(x = c(0,1), size = 1, replace = T, prob = c(0.5,0.5))
                SV_2[m, (n*2):(n*2+1)] <- c(whichhom,whichhom)
              }
            }
          }
        }
        SV_2 <- SV_2[rowSums(SV_2[, 2:n_haps]) > 0, ]
        if(ncol(SV_2)>1 & nrow(SV_2)>1){
          positions_snps <- rbind(positions_snps, cbind(SV_2$V1+pos_counter, rowSums(SV_2[, 2:n_haps]), rowSums(SV_2[, subsample]),"multi-copy"))
          positions_snps_wd <- rbind(positions_snps_wd, cbind(SV_2$V1+pos_counter_wd, rowSums(SV_2[, 2:n_haps]), rowSums(SV_2[, subsample]), "multi-copy"))
        }
      }
      for(j in 2:n_haps){
        new_seq<- seq[pos_counter:(pos_counter+l_SV-1)]
        if(ncol(SV_2)>1){
          derived_pos<-SV_2$V1[SV_2[j]==1]+1
          new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
        }
        fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
      }
      pos_counter <- pos_counter +l_SV
      pos_counter_wd <- pos_counter_wd + l_SV
    }
    print(i)
  }
  SC <- read.table(file = files_SC_sorted[length(files_SC_sorted)], header = F, sep = "")
  l_SC <- l_numbers_SC[length(files_SC_sorted)]
  if(ncol(SC)>1){
    for(m in 1:nrow(SC)){
      for(n in 1:((n_haps-1)/2)){
        if(sum(SC[m, (n*2):(n*2+1)])==1){
          hetorhom <- sample(x = c(1,2), size = 1, replace = T, prob = c((1-Fis), Fis))
          if(hetorhom==2){
            whichhom <- sample(x = c(0,1), size = 1, replace = T, prob = c(0.5,0.5))
            SC[m, (n*2):(n*2+1)] <- c(whichhom,whichhom)
          }
        }
      }
    }
    SC <- SC[rowSums(SC[, 2:n_haps]) > 0, ]
    if(ncol(SC)>1 & nrow(SC)>1){
      positions_snps <- rbind(positions_snps, cbind(SC$V1+pos_counter, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]), "single-copy"))
      positions_snps_wd <- rbind(positions_snps_wd, cbind(SC$V1+pos_counter_wd, rowSums(SC[, 2:n_haps]), rowSums(SC[, subsample]), "single-copy"))
    }
  }
  for(j in 2:n_haps){
    new_seq <- seq[pos_counter:(pos_counter+l_SC-1)]
    if(ncol(SC)>1){
      derived_pos <- SC$V1[SC[j]==1]+1
      new_seq[derived_pos] <- seq_derived[derived_pos+pos_counter-1]
    }
    fastas[[(j-1)]]  <- c(fastas[[(j-1)]],new_seq)
  }

  for(i in 1:length(fastas)){
    fastas[[i]] <- fastas[[i]][!is.na(fastas[[i]])]
  }

  positions_snps <- as.data.frame(positions_snps)
  positions_snps$V1 <- as.numeric(positions_snps$V1)
  positions_snps <- positions_snps[order(positions_snps$V1),]
  # positions_snps <- unique(positions_snps)

  write.table(positions_snps, file = paste0(output_dir, "/Position_files_Fis0.9","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  # Function to wrap sequence into lines of 60 characters

  positions_snps_wd <- as.data.frame(positions_snps_wd)
  positions_snps_wd$V1 <- as.numeric(positions_snps_wd$V1)
  positions_snps_wd <- positions_snps_wd[order(positions_snps_wd$V1),]
  # positions_snps_wd <- unique(positions_snps_wd)

  write.table(positions_snps_wd, file = paste0(output_dir, "/Position_wd_files_Fis0.9","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(subsample-1, file = paste0(output_dir, "/Subsampled_haplotypes_Fis0.9","_Rep",k, ".txt"), sep = "\t", quote = F, row.names = F, col.names = F)

  # Loop over all fastas
  i <-1
  for (i in seq_along(fastas)) {
    header <- fastas[[i]][1]
    sequence <- fastas[[i]][-1]

    # Wrap the sequence
    wrapped_seq <- wrap_sequence(sequence)

    # Define output filename (e.g., Sample1.fasta)
    outfile <- file.path(output_dir, paste0("Fis0.9","_Rep",k,"_Sample", ceiling(i/2), "_hap", (2-((i)%%2)) , ".fasta"))

    # Write to file
    writeLines(c(header, wrapped_seq), con = outfile)
  }
}

