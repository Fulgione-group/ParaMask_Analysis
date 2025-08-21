


sim_stats <- c()
nSNP<-0
i <- 1
i<-2
i<-3
for(i in 1:3){
  #read rCNV results
  rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/HW_prop/rep",as.character(i),"/cnv.txt"))
  nrow(rCNV)
  #only regard first 1 Mbp 
  rCNV <- rCNV[rCNV$POS<=1000000,]
  #simulated copy number info 
  SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/HW_prop/rep",as.character(i),"/Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
  nrow(SVtab)
  #for N=100, 5X segregating sites are equal to ParaMask, rCNV does NOT output unclassified SNPs
  Sim<-read.table(file = paste("$PATH_TO_PARAMASK_OUTPUT/HW_10PercentParalogs/ParaMask_runs/Sim0.1HW_rep",as.character(i),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  nrow(Sim)
  SVtab <- SVtab[SVtab$V2 %in% Sim$Position,]
  nSNP<-nSNP+ nrow(SVtab)/3

  colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

  rCNV$POS<- as.numeric(rCNV$POS)
  SVtab$POS<- as.numeric(SVtab$POS)

  SVtab$rCNV_call <- "U"
  SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat
  a<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  b <- sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  c <- sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  d<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  e<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  f<-sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  g<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)

  SVtab<- SVtab[SVtab$rCNV_call!="U",]
  h<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  i<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  j<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  k<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  l<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)
  sim_stats <- rbind(sim_stats, c(a,b,c,d,e,f,g,h,i,j,k,l))
}
colnames(sim_stats) <- c("Correct_SC", "Incorrect_SC", "Uncertain_SC","Correct_SV", "Incorrect_SV", "Uncertain_SV", "Overall_correct_with_uncertain", "Correct_SC_no_uncertain", "Incorrect_SC_no_uncertain", "Correct_SV_no_uncertain", "Incorrect_SV_no_uncertain", "Overall_correct_no_uncertain")
sim_stats <- rbind(sim_stats, colMeans(sim_stats))
nSNP

nSNP<-0

for(i in 1:3){
  rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(i),"/cnv.txt"))
  rCNV <- rCNV[rCNV$POS<1000000,]
  SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(i),"/Selfing_Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
  nrow(SVtab)
  Sim<-read.table(file = paste("$PATH_TO_PARAMASK_OUTPUT/Fis0.9_10PercentParalogs/ParaMask_runs/Sim0.1Fis0.9_rep",as.character(i),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  SVtab <- SVtab[SVtab$V2 %in% Sim$Position,]
  nSNP<-nSNP+ nrow(SVtab)/3

  colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

  rCNV$POS<- as.numeric(rCNV$POS)
  SVtab$POS<- as.numeric(SVtab$POS)

  SVtab$rCNV_call <- "U"
  SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat
  a<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  b <- sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  c <- sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  d<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  e<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  f<-sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  g<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)

  SVtab<- SVtab[SVtab$rCNV_call!="U",]
  h<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  i<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  j<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  k<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  l<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)
  sim_stats <- rbind(sim_stats, c(a,b,c,d,e,f,g,h,i,j,k,l))
}
sim_stats <- rbind(sim_stats, colMeans(sim_stats[5:7,]))
rownames(sim_stats) <- c(paste0("HW_0.1_rep", 1:3), "HW_average", paste0("Selfing_0.1_rep", 1:3), "Selfing_average")
nSNP

write.table(sim_stats, file = "$PATH_TO_rCNV_OUTPUT/sim_results_rCNV_final.txt", sep = "\t", col.names=T,row.names=T, quote=F)

####cov5 and subsample15
extract_alleles <- function(genotypes) {
  alleles <- unlist(strsplit(genotypes[genotypes != "./."], split = "[/|]"))
  return(alleles)
}


sim_stats <- c()
snpcount<-0
for(i in 1:3){
  rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/cov5_subsample15/HW_prop/sample15_cov5//rep",as.character(i),"/cnv.txt"))
  rCNV <- rCNV[rCNV$POS<=1000000,]
  SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/HW_prop/rep",as.character(i),"/Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
  #GTtab contains the segregating sites after downsampling and reducing coverage
  GTtab<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/cov5_subsample15/HW_prop/sample15_cov5//rep",as.character(i),"/gt.tab.txt"))


  # Determine if a site is segregating
  GTtab$segregating <- apply(GTtab[, grep("^genotype_", names(GTtab))], 1, function(row) {
    alleles <- extract_alleles(row)
    length(unique(alleles)) > 1
  })
  SVtab <- SVtab[SVtab$V2 %in% GTtab$POS[GTtab$segregating],]
  nrow(SVtab)
  SVtab <- SVtab[SVtab$V2<=1000000,]
  SVtab1 <-SVtab

  snpcount <-snpcount+ nrow(SVtab)/3
  colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

  rCNV$POS<- as.numeric(rCNV$POS)
  SVtab$POS<- as.numeric(SVtab$POS)

  SVtab$rCNV_call <- "U"
  SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat
  a<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  b <- sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  c <- sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  d<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  e<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  f<-sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  g<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)

  SVtab<- SVtab[SVtab$rCNV_call!="U",]
  h<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  i<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  j<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  k<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  l<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)
  sim_stats <- rbind(sim_stats, c(a,b,c,d,e,f,g,h,i,j,k,l))
}
snpcount

colnames(sim_stats) <- c("Correct_SC", "Incorrect_SC", "Uncertain_SC","Correct_SV", "Incorrect_SV", "Uncertain_SV", "Overall_correct_with_uncertain", "Correct_SC_no_uncertain", "Incorrect_SC_no_uncertain", "Correct_SV_no_uncertain", "Incorrect_SV_no_uncertain", "Overall_correct_no_uncertain")
sim_stats <- rbind(sim_stats, colMeans(sim_stats))

i<-2

nrow(rCNV)
snpcount<-0
for(i in 1:3){
  rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/cov5_subsample15/selfing/sample15_cov5/rep",as.character(i),"/cnv.txt"))
  rCNV <- rCNV[rCNV$POS<=1000000,]
  SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(i),"/Selfing_Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
  GTtab<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/cov5_subsample15/selfing/sample15_cov5//rep",as.character(i),"/gt.tab.txt"))


  # Determine if a site is segregating
  GTtab$segregating <- apply(GTtab[, grep("^genotype_", names(GTtab))], 1, function(row) {
    alleles <- extract_alleles(row)
    length(unique(alleles)) > 1
  })
  SVtab <- SVtab[SVtab$V2 %in% GTtab$POS[GTtab$segregating],]
  #take only the first 1 Mbp
  SVtab <- SVtab[SVtab$V2<=1000000,]
  snpcount <-snpcount+ nrow(SVtab)/3

  colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

  rCNV$POS<- as.numeric(rCNV$POS)
  SVtab$POS<- as.numeric(SVtab$POS)

  SVtab$rCNV_call <- "U"
  SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat
  a<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  b <- sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  c <- sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  d<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  e<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  f<-sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  g<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)

  SVtab<- SVtab[SVtab$rCNV_call!="U",]
  h<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  i<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
  j<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  k<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
  l<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)
  sim_stats <- rbind(sim_stats, c(a,b,c,d,e,f,g,h,i,j,k,l))
}
sim_stats <- rbind(sim_stats, colMeans(sim_stats[5:7,]))
rownames(sim_stats) <- c(paste0("HW_0.1_rep", 1:3), "HW_average", paste0("Selfing_0.1_rep", 1:3), "Selfing_average")



write.table(sim_stats, file = "$PATH_TO_rCNV_OUTPUT/sim_results_rCNV_final_15N_5X.txt", sep = "\t", col.names=T,row.names=T, quote=F)









#######
rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/rep1_FIS0.9/cnv.txt"))
rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/rep1_FIS/cnv.txt"))

SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(1),"/Selfing_Simulations_SV_SC_rep",as.character(1),"_mat_cpos.biallelic_true.txt"))
colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

rCNV$POS<- as.numeric(rCNV$POS)
SVtab$POS<- as.numeric(SVtab$POS)

SVtab$rCNV_call <- "U"
SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat
a<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
b <- sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
c <- sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
d<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
e<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
f<-sum(SVtab$rCNV_call=="U" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
g<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)

SVtab<- SVtab[SVtab$rCNV_call!="U",]
h<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
i<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SC")/sum(SVtab$true_cn=="SC")
j<-sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
k<-sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SV")/sum(SVtab$true_cn=="SV")
l<-(sum(SVtab$rCNV_call=="cnv" & SVtab$true_cn=="SV") + sum(SVtab$rCNV_call=="non-cnv" & SVtab$true_cn=="SC"))/nrow(SVtab)
sim_stats <- rbind(sim_stats, c(a,b,c,d,e,f,g,h,i,j,k,l))



##plot a single rCNV run
rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/HW_prop/rep",as.character(1),"/cnv.txt"))
SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/HW_prop/rep",as.character(1),"/Simulations_SV_SC_rep",as.character(1),"_mat_cpos.biallelic_true.txt"))
colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")

rCNV$POS<- as.numeric(rCNV$POS)
SVtab$POS<- as.numeric(SVtab$POS)

SVtab$rCNV_call <- "U"
SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat


rCNV$alt_freq<- rCNV$propHet/2 + rCNV$propHomAlt
ncol(rCNV)
colnames(rCNV)
rCNV$MAF <- apply(X = rCNV,MARGIN = 1,FUN = function(x){min(as.numeric(x[12]), (1-as.numeric(x[12])))})
ggplot(data = rCNV)+
  geom_point(aes(x = MAF, y = propHet, color = dup.stat))+
  scale_x_continuous(lim = c(-0.025,0.525), breaks = c(0,0.5), expand=c(0,0))+
  scale_y_continuous(lim = c(-0.05,1.05), breaks = c(0,0.5,1), expand=c(0,0))+
  # scale_linetype_manual(name="SNP", labels=c("single-copy SNP","multicopy SNP"),values=c("solid","twodash"))+
  # scale_color_gradient(name="SNP",limits=c(0,1),breaks=c(1,0), labels=c("single-copy SNP","multicopy SNP"), "",high="blue3", low="brown2")+
  labs(x= "minor allele frequency (maf)", y = "heterozygote frequency")+
  coord_fixed(ratio = 0.5)+
  geom_segment(x=-0.025, y=-0.003, xend=-0.025, yend=1.003, size=1.2)+
  geom_segment(x=-0.0013, y=-0.05, xend=0.5013, yend=-0.05, size=1.2)+
  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=14) , axis.ticks=element_line(size=1),axis.ticks.length = unit(0.25,"cm"), axis.title = element_text(size=16),legend.text = element_text(size=12), legend.position = c(0.2, 0.6),legend.title=element_blank(), legend.key = element_rect(color="white", fill=NA), legend.key.size = unit(1.5,"cm") ,panel.background = element_rect(fill=NA, color="white"))+
  guides(color = guide_legend(override.aes = list(size = 4)))


nrow(rCNV)
rCNV<-read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(1),"/cnv.txt"))
SVtab <- read.table(file = paste0("$PATH_TO_rCNV_OUTPUT/selfing/rep",as.character(1),"/Selfing_Simulations_SV_SC_rep",as.character(1),"_mat_cpos.biallelic_true.txt"))
colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")
nrow(rCNV)
rCNV$POS<- as.numeric(rCNV$POS)
SVtab$POS<- as.numeric(SVtab$POS)

SVtab$rCNV_call <- "U"
SVtab$rCNV_call[SVtab$POS %in% rCNV$POS] <-rCNV$dup.stat

rCNV$alt_freq<- rCNV$propHet/2 + rCNV$propHomAlt
ncol(rCNV)
colnames(rCNV)
rCNV$MAF <- apply(X = rCNV,MARGIN = 1,FUN = function(x){min(as.numeric(x[12]), (1-as.numeric(x[12])))})
ggplot(data = rCNV)+
  geom_point(aes(x = MAF, y = propHet, color = dup.stat))+
  scale_x_continuous(lim = c(-0.025,0.525), breaks = c(0,0.5), expand=c(0,0))+
  scale_y_continuous(lim = c(-0.05,1.05), breaks = c(0,0.5,1), expand=c(0,0))+
  # scale_linetype_manual(name="SNP", labels=c("single-copy SNP","multicopy SNP"),values=c("solid","twodash"))+
  # scale_color_gradient(name="SNP",limits=c(0,1),breaks=c(1,0), labels=c("single-copy SNP","multicopy SNP"), "",high="blue3", low="brown2")+
  labs(x= "minor allele frequency (maf)", y = "heterozygote frequency")+
  coord_fixed(ratio = 0.5)+
  geom_segment(x=-0.025, y=-0.003, xend=-0.025, yend=1.003, size=1.2)+
  geom_segment(x=-0.0013, y=-0.05, xend=0.5013, yend=-0.05, size=1.2)+
  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=14) , axis.ticks=element_line(size=1),axis.ticks.length = unit(0.25,"cm"), axis.title = element_text(size=16),legend.text = element_text(size=12), legend.position = c(0.2, 0.6),legend.title=element_blank(), legend.key = element_rect(color="white", fill=NA), legend.key.size = unit(1.5,"cm") ,panel.background = element_rect(fill=NA, color="white"))+
  guides(color = guide_legend(override.aes = list(size = 4)))
