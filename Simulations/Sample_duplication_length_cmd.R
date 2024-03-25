#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##default values
verbose <- FALSE
outpath <- getwd()
rate <- 1/1000
ID <-""

outpath <- "/home/btjeng/Data/SeDuS/Sim_rep_selfing"


simlen <- 1000000
proportion <- 0.05
proportion <- 0.5
proportion <- 0.3
proportion <- 0.2
proportion <- 0.4
proportion <- 0
proportion <- 0.15
proportion <- 0.25
proportion <- 0.35
proportion <- 0.45
proportion <- 0.1

ID<-"rep1"
ID<-"rep2"
ID<-"rep3"

for(i in 1:length(args)){
  print(args[i])
}

for(i in 1:length(args)){
  if (args[i]=="--simlen" | args[i]=="-sl" ){
    simlen <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--prop" | args[i]=="-p" ){
    proportion <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--rate" | args[i]=="-r" ){
    proportion <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--outpath" | args[i]=="-o" ){
    outpath <- as.character(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--ID" | args[i]=="-id" ){
    ID <- as.character(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--verbose" | args[i]=="-v" ){
    verbose <- TRUE
  } 
}
##
if(proportion > 0){
  seqpl <- proportion*simlen
  
  suml <- 0
  SV_length_sample<- c()
  while(suml < seqpl){
    SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=rate)))
    suml<- sum(SV_length_sample)
  }
  SV_length_sample
  is_truncate<-sample(x = c(TRUE,FALSE), size = 1)
  if(is_truncate){
    tmp<-sample(c(1:length(SV_length_sample)), size = 1)
    suml <- suml - SV_length_sample[tmp]
    SV_length_sample <- SV_length_sample[-tmp]
  }
  
  seqscl<-simlen-suml
  
  
  sumscl<-0
  
  # SC_length_sample<-rgeom(1, prob = length(SV_length_sample)/seqscl)
  
  #seqscl
  SC_length_sample<- c()
  
  seqscl_rest <- seqscl
  
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)+1)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[1]
  for(i in 1:length(SV_length_sample)){
    SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+1)/seqscl_rest))
    sumscl<- sum(SC_length_sample)
    seqscl_rest<- seqscl_rest-SC_length_sample[(i+1)]
  }
  tmpEND<-SC_length_sample[length(SC_length_sample)]
  while(is.na(tmpEND)){
    SC_length_sample<- c()
    seqscl_rest <- seqscl
    
    SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)+1)/seqscl_rest))
    sumscl<- sum(SC_length_sample)
    seqscl_rest<- seqscl_rest-SC_length_sample[1]
    for(i in 1:length(SV_length_sample)){
      SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+1)/seqscl_rest))
      sumscl<- sum(SC_length_sample)
      seqscl_rest<- seqscl_rest-SC_length_sample[(i+1)]
    }
    tmpEND<-SC_length_sample[length(SC_length_sample)]
  }
  # sum(SC_length_sample)
  # SC_length_sample
  
  # there was an error here
  SC_length_sample[length(SC_length_sample)]<- SC_length_sample[length(SC_length_sample)] + seqscl_rest
  # sum(SC_length_sample)+ sum(SV_length_sample)
  
  # length(SC_length_sample)
  # length(SV_length_sample)
  # head(SC_length_sample)
  
  ##convert to positions
  SVtab <- c()
  pos<- SC_length_sample[1]
  SVtab<- c(1, pos, "sc")
  i<-1
  for(i in 1:length(SV_length_sample)){
    pos2<-pos + SV_length_sample[i]
    SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
    pos <- pos2 
    pos2<- pos + SC_length_sample[i+1]
    SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
    pos <- pos2 
  }
  SVtab <- data.frame(SVtab,stringsAsFactors = F)
  colnames(SVtab)<- c("start", "end", "type")
  
  #not needed anymore
  #SVtab$end[nrow(SVtab)] <- simlen
  head(SVtab)
  tail(SVtab)

  
  
  write.table(x = SVtab, file = paste(outpath, "/SVtab_", as.character(proportion), "_",ID,".txt", sep = ""), sep = "\t", row.names = F, quote = F, col.names = T)
  write.table(x = SC_length_sample, file = paste(outpath, "/SC_length_", as.character(proportion), "_",ID,".txt", sep = ""), sep = "\t", row.names = F, quote = F, col.names = F)
  write.table(x = SV_length_sample, file = paste(outpath, "/SV_length_", as.character(proportion), "_",ID,".txt", sep = ""), row.names = F, quote = F, col.names = F)
} else if (proportion==0){

  
  SC_length_sample<- c()
  seqscl_rest <- simlen
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (200+1)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[1]
  for(i in 1:200){
    SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (200-i+1)/seqscl_rest))
    sumscl<- sum(SC_length_sample)
    seqscl_rest<- seqscl_rest-SC_length_sample[(i+1)]
  }
  # sum(SC_length_sample)
  # SC_length_sample
  
  # there was an error here
  SC_length_sample[length(SC_length_sample)]<- SC_length_sample[length(SC_length_sample)] + seqscl_rest
  SVtab <- c()
  pos<- SC_length_sample[1]
  SVtab<- c(1, pos, "sc")
  i<-1
  for(i in 1:(length(SC_length_sample)-1)){
    pos2<- pos + SC_length_sample[i+1]
    SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
    pos <- pos2 
  }
  SVtab <- data.frame(SVtab,stringsAsFactors = F)
  colnames(SVtab)<- c("start", "end", "type")
  write.table(x = SVtab, file = paste(outpath, "/SVtab_", as.character(proportion), "_",ID,".txt", sep = ""), sep = "\t", row.names = F, quote = F, col.names = T)
  write.table(x = SC_length_sample, file = paste(outpath, "/SC_length_", as.character(proportion), "_",ID,".txt", sep = ""), sep = "\t", row.names = F, quote = F, col.names = F)
}
  
