library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpattern)


#######calculate Overlap of ParaMask with simulations across replicates
result_tab<-c()
reps <- 1:3
k<-1
for(k in reps){
  print(k)
  Sim<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_10PercentParalogs/ParaMask_runs/Sim0.1HW_rep",as.character(k),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  SVtab<- read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_10PercentParalogs/SVtab_0.1_rep",as.character(k) ,".txt", sep=""), header = T, sep = "\t")
  SVtab$cluster <- 0
  SVtab$cluster[SVtab$type=="par"] <-1:length(SVtab$cluster[SVtab$type=="par"])
  
  head(Sim)
  head(SVtab)
  
  Sim$simulated_type <- NA
  Sim$sim_cluster<-0
  ccounter<-0
  prev <- 0
  i<-1
  c<-1
  for(i in 1:nrow(Sim)){
    snp_pos <- Sim$Position[i]
    while(snp_pos <= SVtab$end[nrow(SVtab)]){
      if(snp_pos >= SVtab$start[c] & snp_pos <= SVtab$end[c]){
        Sim$simulated_type[i] <- SVtab$type[c]
        if(Sim$simulated_type[i]!=prev){
          ccounter<- ccounter + 1
        }
        if(Sim$simulated_type[i]=="par"){
          Sim$sim_cluster[i] <-ccounter
          ccounter<-ccounter
        }
        prev<- Sim$simulated_type[i]
        break
      }else{
        c <- c + 1
      }
    }
  }
  Sim$sim_cluster<- Sim$sim_cluster/2
  uc <-unique(Sim$sim_cluster)[-1]
  
  tail(Sim)
  per_par_recover<-c()
  for(i in 1:length(uc)){
    tmp <- Sim[Sim$sim_cluster==uc[i],]
    per_par_recover <-rbind( per_par_recover, c(sum(tmp$finalClass==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_recover<- as.data.frame(per_par_recover)
  colnames(per_par_recover)<- c("recovery", "nSNPs")
  per_par_recover$weight <- per_par_recover$nSNPs/sum(per_par_recover$nSNPs)
  sum(per_par_recover$recovery* per_par_recover$weight)
  median(per_par_recover$recovery)
  mean(per_par_recover$recovery)
  
  
  
  
  #sc overlap
  total_sc_snps <- sum(Sim$simulated_type=="sc")
  total_par_snps <- sum(Sim$simulated_type=="par")
  total_sc_overlap<- sum(Sim$simulated_type=="sc" & Sim$finalClass==0 )
  
  total_wrong_sc<- sum(Sim$simulated_type=="sc" & Sim$finalClass==1 )

  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  
      
  total_sc_overlap /total_sc_snps
  
  #SV overlap
  total_par_overlap <- sum(Sim$simulated_type=="par" & Sim$finalClass==1 )
  
  total_wrong_par <- sum(Sim$simulated_type=="par" & Sim$finalClass==0 )
  
  
  total_par_overlap/total_par_snps
  
  #check how many clusters are picked up
  SVtab$nSNPs <-0
  SVtab$nParalogSNPs <-0
  
  i<-2
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      SVtab$nSNPs[i] <- nrow(Sim[Sim$Position>=SVtab$start[i] &Sim$Position<=SVtab$end[i],])
      SVtab$nParalogSNPs[i] <-nrow(Sim[Sim$Position>=SVtab$start[i] & Sim$Position<=SVtab$end[i] & Sim$finalClass==1,])
    }
  }
  
  
  sum(SVtab$type=="par")
  sum(SVtab$nSNPs>0)
  sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)
  
  #
  
  SVtab_infer<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_10PercentParalogs/ParaMask_runs/Sim0.1HW_rep",as.character(k),"EMresults.finalClass.bed", sep = ""), header = T, sep = "\t")
  
  
  overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000), cluster=rep(0, 1000000))
  head(overlap)
  for(i in 1:nrow(SVtab_infer)){
    if(SVtab_infer$type[i]==1){
      overlap$inferred[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- 1
    }
  }
  
  ccounter <-1
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <-1
      overlap$cluster[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]]<- ccounter
      ccounter<-ccounter + 1
    }
  }
  
  
  
  head(overlap)
  overlap$overlap<-0
  overlap$overlap[overlap$real==0 & overlap$inferred==0]<- 1
  overlap$overlap[overlap$real==1 & overlap$inferred==1]<- 2
  overlap$overlap[overlap$real==1 & overlap$inferred==0]<- 3
  
  uc <-  unique(overlap$cluster)[-1]
  per_par_agreement <-c()
  for(i in 1:length(uc)){
    tmp <- overlap[overlap$cluster==uc[i],]
    per_par_agreement <-rbind(per_par_agreement, c(sum(tmp$inferred==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_agreement<- as.data.frame(per_par_agreement)
  colnames(per_par_agreement) <-c("agreement", "nSNPs")
  per_par_agreement$weight <- per_par_agreement$nSNPs/sum(per_par_agreement$nSNPs)
  sum(per_par_agreement$agreement*per_par_agreement$weight)
  mean(per_par_agreement$agreement)
  median(per_par_agreement$agreement)
  
  sc_agreement <- sum(overlap$overlap==1)/sum(overlap$real==0)
  sc_discrepency <- sum(overlap$overlap==0)/sum(overlap$real==0)
  
  par_agreement <- sum(overlap$overlap==2)/sum(overlap$real==1)
  par_discrepency <- sum(overlap$overlap==3)/sum(overlap$real==1)
  
  overall_agreement <- sum(overlap$overlap==1 | overlap$overlap==2)/nrow(overlap)
  
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_par_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="par")
  
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="sc")
  total_wrong_sc_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="sc") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="par") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  total_wrong_SC_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="sc")
  
  
  result_tab <- rbind(result_tab,c(0.1, 0, total_sc_snps, total_sc_overlap /total_sc_snps, total_wrong_sc/total_sc_snps, total_par_snps,total_par_overlap/total_par_snps, total_wrong_par/total_par_snps,(total_sc_overlap + total_par_overlap)/(total_sc_snps+total_par_snps),total_correct_par_EM_snps/total_par_snps,total_wrong_par_EM_snps/total_par_snps, total_correct_sc_EM_snps/total_sc_snps, total_wrong_sc_EM_snps/total_sc_snps,total_correct_EM_snps/(total_sc_snps +total_par_snps), total_wrong_EM_snps/(total_sc_snps +total_par_snps), total_correct_par_AR_snps/total_par_snps, total_wrong_SC_AR_snps/total_sc_snps,sc_agreement, sc_discrepency, par_agreement, par_discrepency, overall_agreement, sum(SVtab$type=="par"),sum(SVtab$nSNPs>0),sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)))
  
}


result_tab


reps <- 1:3
k<-1
for(k in reps){
  print(k)
  Sim<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_50PercentParalogs/ParaMask_runs/Sim0.5HW_rep",as.character(k),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  SVtab<- read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_50PercentParalogs/SVtab_0.5_rep",as.character(k) ,".txt", sep=""), header = T, sep = "\t")
  
  head(Sim)
  head(SVtab)
  
  Sim$simulated_type <- NA
  Sim$sim_cluster<-0
  ccounter<-0
  prev <- 0
  i<-1
  c<-1
  for(i in 1:nrow(Sim)){
    snp_pos <- Sim$Position[i]
    while(snp_pos <= SVtab$end[nrow(SVtab)]){
      if(snp_pos >= SVtab$start[c] & snp_pos <= SVtab$end[c]){
        Sim$simulated_type[i] <- SVtab$type[c]
        if(Sim$simulated_type[i]!=prev){
          ccounter<- ccounter + 1
        }
        if(Sim$simulated_type[i]=="par"){
          Sim$sim_cluster[i] <-ccounter
          ccounter<-ccounter
        }
        prev<- Sim$simulated_type[i]
        break
      }else{
        c <- c + 1
      }
    }
  }
  Sim$sim_cluster<- Sim$sim_cluster/2
  uc <-unique(Sim$sim_cluster)[-1]
  
  per_par_recover<-c()
  for(i in 1:length(uc)){
    tmp <- Sim[Sim$sim_cluster==uc[i],]
    per_par_recover <-rbind( per_par_recover, c(sum(tmp$finalClass==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_recover<- as.data.frame(per_par_recover)
  colnames(per_par_recover)<- c("recovery", "nSNPs")
  per_par_recover$weight <- per_par_recover$nSNPs/sum(per_par_recover$nSNPs)
  sum(per_par_recover$recovery* per_par_recover$weight)
  median(per_par_recover$recovery)
  mean(per_par_recover$recovery)
  
  total_sc_snps<-sum(Sim$simulated_type=="sc")
  total_par_snps<-sum(Sim$simulated_type=="par")
  
  #sc overlap
  total_sc_overlap<- sum(Sim$simulated_type=="sc" & Sim$finalClass==0 )
  
  total_wrong_sc<- sum(Sim$simulated_type=="sc" & Sim$finalClass==1 )
  
  
  total_sc_overlap/total_sc_snps
  #SV overlap
  total_par_overlap <- sum(Sim$simulated_type=="par" & Sim$finalClass==1 )
  
  total_wrong_par <- sum(Sim$simulated_type=="par" & Sim$finalClass==0 )
  total_par_overlap/total_par_snps
  #check how many clusters are picked up
  SVtab$nSNPs <-0
  SVtab$nParalogSNPs <-0
  
  i<-2
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      SVtab$nSNPs[i] <- nrow(Sim[Sim$Position>=SVtab$start[i] &Sim$Position<=SVtab$end[i],])
      SVtab$nParalogSNPs[i] <-nrow(Sim[Sim$Position>=SVtab$start[i] & Sim$Position<=SVtab$end[i] & Sim$finalClass==1,])
    }
  }
  
  
  sum(SVtab$type=="par")
  sum(SVtab$nSNPs>0)
  sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)
  
  #
  
  
  SVtab_infer<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/HW_50PercentParalogs/ParaMask_runs/Sim0.5HW_rep",as.character(k),"EMresults.finalClass.bed", sep = ""), header = T, sep = "\t")
  
  
  overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000), cluster=rep(0, 1000000))
  head(overlap)
  for(i in 1:nrow(SVtab_infer)){
    if(SVtab_infer$type[i]==1){
      overlap$inferred[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- 1
    }
  }
  
  ccounter <-1
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <-1
      overlap$cluster[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]]<- ccounter
      ccounter<-ccounter + 1
    }
  }
  

  
  head(overlap)
  overlap$overlap<-0
  overlap$overlap[overlap$real==0 & overlap$inferred==0]<- 1
  overlap$overlap[overlap$real==1 & overlap$inferred==1]<- 2
  overlap$overlap[overlap$real==1 & overlap$inferred==0]<- 3

  uc <-  unique(overlap$cluster)[-1]
  per_par_agreement <-c()
  for(i in 1:length(uc)){
    tmp <- overlap[overlap$cluster==uc[i],]
    per_par_agreement <-rbind(per_par_agreement, c(sum(tmp$inferred==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_agreement<- as.data.frame(per_par_agreement)
  colnames(per_par_agreement) <-c("agreement", "nSNPs")
  per_par_agreement$weight <- per_par_agreement$nSNPs/sum(per_par_agreement$nSNPs)
  sum(per_par_agreement$agreement*per_par_agreement$weight)
  mean(per_par_agreement$agreement)
  median(per_par_agreement$agreement)
  
  sc_agreement <- sum(overlap$overlap==1)/sum(overlap$real==0)
  sc_discrepency <- sum(overlap$overlap==0)/sum(overlap$real==0)
  
  par_agreement <- sum(overlap$overlap==2)/sum(overlap$real==1)
  par_discrepency <- sum(overlap$overlap==3)/sum(overlap$real==1)
  
  ##
  
  
  ##
  
  overall_agreement <- sum(overlap$overlap==1 | overlap$overlap==2)/nrow(overlap)
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_par_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="par")
  
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="sc")
  total_wrong_sc_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="sc") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="par") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  total_wrong_SC_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="sc")
  
  
  result_tab <- rbind(result_tab,c(0.5, 0, total_sc_snps, total_sc_overlap /total_sc_snps, total_wrong_sc/total_sc_snps, total_par_snps,total_par_overlap/total_par_snps, total_wrong_par/total_par_snps,(total_sc_overlap + total_par_overlap)/(total_sc_snps+total_par_snps),total_correct_par_EM_snps/total_par_snps,total_wrong_par_EM_snps/total_par_snps, total_correct_sc_EM_snps/total_sc_snps, total_wrong_sc_EM_snps/total_sc_snps,total_correct_EM_snps/(total_sc_snps +total_par_snps), total_wrong_EM_snps/(total_sc_snps +total_par_snps), total_correct_par_AR_snps/total_par_snps, total_wrong_SC_AR_snps/total_sc_snps,sc_agreement, sc_discrepency, par_agreement, par_discrepency, overall_agreement, sum(SVtab$type=="par"),sum(SVtab$nSNPs>0),sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)))
  
}


reps <- 1:3
k<-3
for(k in reps){
  print(k)
  Sim<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_10PercentParalogs/ParaMask_runs/Sim0.1Fis0.9_rep",as.character(k),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  SVtab<- read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_10PercentParalogs/SVtab_0.1_rep",as.character(k) ,".txt", sep=""), header = T, sep = "\t")
  SVtab$cluster <- 0
  SVtab$cluster[SVtab$type=="par"] <-1:length(SVtab$cluster[SVtab$type=="par"])
  
  head(Sim)
  head(SVtab)
  
  Sim$simulated_type <- NA
  Sim$sim_cluster<-0
  ccounter<-0
  prev <- 0
  i<-1
  c<-1
  for(i in 1:nrow(Sim)){
    snp_pos <- Sim$Position[i]
    while(snp_pos <= SVtab$end[nrow(SVtab)]){
      if(snp_pos >= SVtab$start[c] & snp_pos <= SVtab$end[c]){
        Sim$simulated_type[i] <- SVtab$type[c]
        if(Sim$simulated_type[i]!=prev){
          ccounter<- ccounter + 1
        }
        if(Sim$simulated_type[i]=="par"){
          Sim$sim_cluster[i] <-ccounter
          ccounter<-ccounter
        }
        prev<- Sim$simulated_type[i]
        break
      }else{
        c <- c + 1
      }
    }
  }
  Sim$sim_cluster<- Sim$sim_cluster/2
  uc <-unique(Sim$sim_cluster)[-1]
  
  per_par_recover<-c()
  for(i in 1:length(uc)){
    tmp <- Sim[Sim$sim_cluster==uc[i],]
    per_par_recover <-rbind( per_par_recover, c(sum(tmp$finalClass==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_recover<- as.data.frame(per_par_recover)
  colnames(per_par_recover)<- c("recovery", "nSNPs")
  per_par_recover$weight <- per_par_recover$nSNPs/sum(per_par_recover$nSNPs)
  sum(per_par_recover$recovery* per_par_recover$weight)
  median(per_par_recover$recovery)
  mean(per_par_recover$recovery)
  
  total_sc_snps <- sum(Sim$simulated_type=="sc")
  total_par_snps <- sum(Sim$simulated_type=="par")
  
  
  
  #sc overlap
  total_sc_overlap<- sum(Sim$simulated_type=="sc" & Sim$finalClass==0 )
  
  total_wrong_sc<- sum(Sim$simulated_type=="sc" & Sim$finalClass==1 )
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  
  
  total_sc_overlap /total_sc_snps
  
  #SV overlap
  total_par_overlap <- sum(Sim$simulated_type=="par" & Sim$finalClass==1 )
  
  total_wrong_par <- sum(Sim$simulated_type=="par" & Sim$finalClass==0 )
  
  
  total_par_overlap/total_par_snps
  
  #check how many clusters are picked up
  SVtab$nSNPs <-0
  SVtab$nParalogSNPs <-0
  
  i<-2
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      SVtab$nSNPs[i] <- nrow(Sim[Sim$Position>=SVtab$start[i] &Sim$Position<=SVtab$end[i],])
      SVtab$nParalogSNPs[i] <-nrow(Sim[Sim$Position>=SVtab$start[i] & Sim$Position<=SVtab$end[i] & Sim$finalClass==1,])
    }
  }
  
  
  sum(SVtab$type=="par")
  sum(SVtab$nSNPs>0)
  sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)
  
  #
  
  SVtab_infer<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_10PercentParalogs/ParaMask_runs/Sim0.1Fis0.9_rep",as.character(k),"EMresults.finalClass.bed", sep = ""), header = T, sep = "\t")
  
  
  overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000), cluster=rep(0, 1000000))
  head(overlap)
  for(i in 1:nrow(SVtab_infer)){
    if(SVtab_infer$type[i]==1){
      overlap$inferred[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- 1
    }
  }
  
  ccounter <-1
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <-1
      overlap$cluster[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]]<- ccounter
      ccounter<-ccounter + 1
    }
  }
  
  
  
  head(overlap)
  overlap$overlap<-0
  overlap$overlap[overlap$real==0 & overlap$inferred==0]<- 1
  overlap$overlap[overlap$real==1 & overlap$inferred==1]<- 2
  overlap$overlap[overlap$real==1 & overlap$inferred==0]<- 3
  
  uc <-  unique(overlap$cluster)[-1]
  per_par_agreement <-c()
  for(i in 1:length(uc)){
    tmp <- overlap[overlap$cluster==uc[i],]
    per_par_agreement <-rbind(per_par_agreement, c(sum(tmp$inferred==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_agreement<- as.data.frame(per_par_agreement)
  colnames(per_par_agreement) <-c("agreement", "nSNPs")
  per_par_agreement$weight <- per_par_agreement$nSNPs/sum(per_par_agreement$nSNPs)
  sum(per_par_agreement$agreement*per_par_agreement$weight)
  mean(per_par_agreement$agreement)
  median(per_par_agreement$agreement)
  
  sc_agreement <- sum(overlap$overlap==1)/sum(overlap$real==0)
  sc_discrepency <- sum(overlap$overlap==0)/sum(overlap$real==0)
  
  par_agreement <- sum(overlap$overlap==2)/sum(overlap$real==1)
  par_discrepency <- sum(overlap$overlap==3)/sum(overlap$real==1)
  
  overall_agreement <- sum(overlap$overlap==1 | overlap$overlap==2)/nrow(overlap)
  
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_par_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="par")
  
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="sc")
  total_wrong_sc_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="sc") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="par") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  total_wrong_SC_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="sc")
  
  
  result_tab <- rbind(result_tab,c(0.1, 0.9, total_sc_snps, total_sc_overlap /total_sc_snps, total_wrong_sc/total_sc_snps, total_par_snps,total_par_overlap/total_par_snps, total_wrong_par/total_par_snps,(total_sc_overlap + total_par_overlap)/(total_sc_snps+total_par_snps),total_correct_par_EM_snps/total_par_snps,total_wrong_par_EM_snps/total_par_snps, total_correct_sc_EM_snps/total_sc_snps, total_wrong_sc_EM_snps/total_sc_snps,total_correct_EM_snps/(total_sc_snps +total_par_snps), total_wrong_EM_snps/(total_sc_snps +total_par_snps), total_correct_par_AR_snps/total_par_snps, total_wrong_SC_AR_snps/total_sc_snps,sc_agreement, sc_discrepency, par_agreement, par_discrepency, overall_agreement, sum(SVtab$type=="par"),sum(SVtab$nSNPs>0),sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)))
  
}

reps <- 1:3
k<-2
for(k in reps){
  print(k)
  Sim<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_50PercentParalogs/ParaMask_runs/Sim0.5Fis0.9_rep",as.character(k),"EMresults.finalClass.het", sep = ""), header = T, sep = "\t")
  SVtab<- read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_50PercentParalogs/SVtab_0.5_rep",as.character(k) ,".txt", sep=""), header = T, sep = "\t")
  SVtab$cluster <- 0
  SVtab$cluster[SVtab$type=="par"] <-1:length(SVtab$cluster[SVtab$type=="par"])
  
  head(Sim)
  head(SVtab)
  
  Sim$simulated_type <- NA
  Sim$sim_cluster<-0
  ccounter<-0
  prev <- 0
  i<-1
  c<-1
  for(i in 1:nrow(Sim)){
    snp_pos <- Sim$Position[i]
    while(snp_pos <= SVtab$end[nrow(SVtab)]){
      if(snp_pos >= SVtab$start[c] & snp_pos <= SVtab$end[c]){
        Sim$simulated_type[i] <- SVtab$type[c]
        if(Sim$simulated_type[i]!=prev){
          ccounter<- ccounter + 1
        }
        if(Sim$simulated_type[i]=="par"){
          Sim$sim_cluster[i] <-ccounter
          ccounter<-ccounter
        }
        prev<- Sim$simulated_type[i]
        break
      }else{
        c <- c + 1
      }
    }
  }
  Sim$sim_cluster<- Sim$sim_cluster/2
  uc <-unique(Sim$sim_cluster)[-1]
  
  per_par_recover<-c()
  for(i in 1:length(uc)){
    tmp <- Sim[Sim$sim_cluster==uc[i],]
    per_par_recover <-rbind( per_par_recover, c(sum(tmp$finalClass==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_recover<- as.data.frame(per_par_recover)
  colnames(per_par_recover)<- c("recovery", "nSNPs")
  per_par_recover$weight <- per_par_recover$nSNPs/sum(per_par_recover$nSNPs)
  sum(per_par_recover$recovery* per_par_recover$weight)
  median(per_par_recover$recovery)
  mean(per_par_recover$recovery)
  
  
  total_sc_snps <- sum(Sim$simulated_type=="sc")
  total_par_snps <- sum(Sim$simulated_type=="par")
  
  
  #sc overlap
  total_sc_overlap<- sum(Sim$simulated_type=="sc" & Sim$finalClass==0 )
  
  total_wrong_sc<- sum(Sim$simulated_type=="sc" & Sim$finalClass==1 )
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  
  
  total_sc_overlap /total_sc_snps
  
  #SV overlap
  total_par_overlap <- sum(Sim$simulated_type=="par" & Sim$finalClass==1 )
  
  total_wrong_par <- sum(Sim$simulated_type=="par" & Sim$finalClass==0 )
  
  
  total_par_overlap/total_par_snps
  
  #check how many clusters are picked up
  SVtab$nSNPs <-0
  SVtab$nParalogSNPs <-0
  
  i<-2
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      SVtab$nSNPs[i] <- nrow(Sim[Sim$Position>=SVtab$start[i] &Sim$Position<=SVtab$end[i],])
      SVtab$nParalogSNPs[i] <-nrow(Sim[Sim$Position>=SVtab$start[i] & Sim$Position<=SVtab$end[i] & Sim$finalClass==1,])
    }
  }
  
  
  sum(SVtab$type=="par")
  sum(SVtab$nSNPs>0)
  sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)
  
  #
  
  SVtab_infer<-read.table(file = paste("/home/btjeng/Data/Paramask/Sim_final/Fis0.9_50PercentParalogs/ParaMask_runs/Sim0.5Fis0.9_rep",as.character(k),"EMresults.finalClass.bed", sep = ""), header = T, sep = "\t")
  
  
  overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000), cluster=rep(0, 1000000))
  head(overlap)
  for(i in 1:nrow(SVtab_infer)){
    if(SVtab_infer$type[i]==1){
      overlap$inferred[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- 1
    }
  }
  
  ccounter <-1
  for(i in 1:nrow(SVtab)){
    if(SVtab$type[i]=="par"){
      overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <-1
      overlap$cluster[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]]<- ccounter
      ccounter<-ccounter + 1
    }
  }
  
  
  
  head(overlap)
  overlap$overlap<-0
  overlap$overlap[overlap$real==0 & overlap$inferred==0]<- 1
  overlap$overlap[overlap$real==1 & overlap$inferred==1]<- 2
  overlap$overlap[overlap$real==1 & overlap$inferred==0]<- 3
  
  uc <-  unique(overlap$cluster)[-1]
  per_par_agreement <-c()
  for(i in 1:length(uc)){
    tmp <- overlap[overlap$cluster==uc[i],]
    per_par_agreement <-rbind(per_par_agreement, c(sum(tmp$inferred==1)/nrow(tmp), nrow(tmp))) 
  }
  per_par_agreement<- as.data.frame(per_par_agreement)
  colnames(per_par_agreement) <-c("agreement", "nSNPs")
  per_par_agreement$weight <- per_par_agreement$nSNPs/sum(per_par_agreement$nSNPs)
  sum(per_par_agreement$agreement*per_par_agreement$weight)
  mean(per_par_agreement$agreement)
  median(per_par_agreement$agreement)
  
  sc_agreement <- sum(overlap$overlap==1)/sum(overlap$real==0)
  sc_discrepency <- sum(overlap$overlap==0)/sum(overlap$real==0)
  
  par_agreement <- sum(overlap$overlap==2)/sum(overlap$real==1)
  par_discrepency <- sum(overlap$overlap==3)/sum(overlap$real==1)
  
  overall_agreement <- sum(overlap$overlap==1 | overlap$overlap==2)/nrow(overlap)
  
  
  total_correct_par_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_par_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="par")
  
  total_correct_sc_EM_snps<-sum(Sim$EM_class==0  & Sim$simulated_type=="sc")
  total_wrong_sc_EM_snps<-sum(Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="sc") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="par")
  total_wrong_EM_snps<-sum((Sim$EM_class==0  & Sim$simulated_type=="par") | Sim$EM_class==2 & Sim$allele.deviation.seed==0 & Sim$simulated_type=="sc")
  
  total_correct_par_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="par")
  total_wrong_SC_AR_snps<-sum(Sim$allele.deviation.seed==1 & Sim$simulated_type=="sc")

  
  result_tab <- rbind(result_tab,c(0.5, 0.9, total_sc_snps, total_sc_overlap /total_sc_snps, total_wrong_sc/total_sc_snps, total_par_snps,total_par_overlap/total_par_snps, total_wrong_par/total_par_snps,(total_sc_overlap + total_par_overlap)/(total_sc_snps+total_par_snps),total_correct_par_EM_snps/total_par_snps,total_wrong_par_EM_snps/total_par_snps, total_correct_sc_EM_snps/total_sc_snps, total_wrong_sc_EM_snps/total_sc_snps,total_correct_EM_snps/(total_sc_snps +total_par_snps), total_wrong_EM_snps/(total_sc_snps +total_par_snps), total_correct_par_AR_snps/total_par_snps, total_wrong_SC_AR_snps/total_sc_snps,sc_agreement, sc_discrepency, par_agreement, par_discrepency, overall_agreement, sum(SVtab$type=="par"),sum(SVtab$nSNPs>0),sum(SVtab$nSNPs>0 & SVtab$nParalogSNPs >0)))
  
}

result_tab <- as.data.frame(result_tab)
colnames(result_tab) <- c("Paralog_prop", "Fis","total_SC_SNPs", "SC_correctlyCall_prop", "SC_notCalled_prop", "Total_Par_SNPs", "Par_correctlyCall_prop", "Par_notCalled_prop",  "Called_prop","Par_CalledEM_prop", "Par_CalledwrongEM_prop",  "SC_CalledEM_prop", "SC_CalledwrongEM_prop",   "Total_correctEM_prop",   "Total_wrongEM_prop" , "Par_CalledAR_prop", "SC_CalledwrongAR_prop", "Agreement_in_SingleCopyRegions", "disagreement_in_SingleCopyRegions",  "Agreement_in_MultiCopyRegions", "disagreement_in_MultiCopyRegions", "Overall_Agreement", "Paralogs", "ParalogsWithSNPs", "ParalogTaggedbyParaMask")

# write.table(x = result_tab,file = "~/Data/Paramask/Sim_final/results_HW_Fis0.9_10P_50P.txt", col.names = T, quote = F, sep="\t")
# result_tab<- read.table(file = "~/Data/Paramask/Sim_final/results_HW_Fis0.9_10P_50P.txt",header = T, sep = "\t")

head(result_tab)
result_tab$Total_correctEM_prop[1:3]
result_tab$Total_wrongEM_prop[1:3]

#total SNPs
mean(result_tab$total_SC_SNPs[1:3]+result_tab$Total_Par_SNPs[1:3])
mean(result_tab$Total_Par_SNPs[1:3])

igc_F10OP<- c(69,58,52)
mean(result_tab$total_SC_SNPs[7:9]+result_tab$Total_Par_SNPs[7:9])
igc_F10OP/result_tab$Total_Par_SNPs[7:9]
mean(igc_F10OP/result_tab$Total_Par_SNPs[7:9])

igc_HW10OP<- c(122,154,113)
igc_HW10OP/result_tab$Total_Par_SNPs[1:3]
mean(igc_HW10OP/result_tab$Total_Par_SNPs[1:3])

mean(result_tab$total_SC_SNPs[1:3])
mean(result_tab$Total_Par_SNPs[1:3])
mean(result_tab$total_SC_SNPs[4:6])
mean(result_tab$Total_Par_SNPs[4:6])
mean(result_tab$total_SC_SNPs[4:6]+result_tab$Total_Par_SNPs[4:6])

igc_HW50P <- c(665,619,631)
igc_HW50P/result_tab$Total_Par_SNPs[4:6]
mean(igc_HW50P/result_tab$Total_Par_SNPs[4:6])

mean(result_tab$total_SC_SNPs[7:9])
mean(result_tab$Total_Par_SNPs[7:9])
mean(result_tab$total_SC_SNPs[10:12])
mean(result_tab$Total_Par_SNPs[10:12])
mean(result_tab$total_SC_SNPs[10:12]+result_tab$Total_Par_SNPs[10:12])

igc_F50P <- c(331, 320, 366)
igc_F50P/result_tab$Total_Par_SNPs[10:12]
mean(igc_F50P/result_tab$Total_Par_SNPs[10:12])
##how many pars without SNP
mean(1-result_tab$ParalogsWithSNPs[1:3]/result_tab$Paralogs[1:3])
mean(1-result_tab$ParalogsWithSNPs[4:6]/result_tab$Paralogs[4:6])
mean(1-result_tab$ParalogsWithSNPs[7:9]/result_tab$Paralogs[7:9])
mean(1-result_tab$ParalogsWithSNPs[10:12]/result_tab$Paralogs[10:12])

##EM

mean(result_tab$SC_CalledEM_prop[1:3])
mean(result_tab$SC_CalledwrongEM_prop[1:3])
mean(1-result_tab$SC_CalledwrongEM_prop[1:3]-result_tab$SC_CalledEM_prop[1:3])


mean(result_tab$Par_CalledEM_prop[1:3])
mean(result_tab$Par_CalledwrongEM_prop[1:3])
mean(1-result_tab$Par_CalledwrongEM_prop[1:3]-result_tab$Par_CalledEM_prop[1:3])

mean(result_tab$Par_CalledAR_prop[1:3])

#


mean(result_tab$SC_CalledEM_prop[4:6])
mean(result_tab$SC_CalledwrongEM_prop[4:6])
mean(1-result_tab$SC_CalledwrongEM_prop[4:6]-result_tab$SC_CalledEM_prop[4:6])

mean(result_tab$Par_CalledEM_prop[4:6])
mean(result_tab$Par_CalledwrongEM_prop[4:6])
mean(1-result_tab$Par_CalledwrongEM_prop[4:6]-result_tab$Par_CalledEM_prop[4:6])

mean(result_tab$Par_CalledAR_prop[4:6])

#
mean(result_tab$SC_CalledEM_prop[7:9])
mean(result_tab$SC_CalledwrongEM_prop[7:9])
mean(1-result_tab$SC_CalledwrongEM_prop[7:9]-result_tab$SC_CalledEM_prop[7:9])

mean(result_tab$Par_CalledEM_prop[7:9])
mean(result_tab$Par_CalledwrongEM_prop[7:9])
mean(1-result_tab$Par_CalledwrongEM_prop[7:9]-result_tab$Par_CalledEM_prop[7:9])


mean(result_tab$Par_CalledAR_prop[1:3]/result_tab$Par_CalledAR_prop[7:9])

#
mean(result_tab$SC_CalledEM_prop[10:12])
mean(result_tab$SC_CalledwrongEM_prop[10:12])
mean(1-result_tab$SC_CalledwrongEM_prop[10:12]-result_tab$SC_CalledEM_prop[10:12])

mean(result_tab$Par_CalledEM_prop[10:12])
mean(result_tab$Par_CalledwrongEM_prop[10:12])
mean(1-result_tab$Par_CalledwrongEM_prop[10:12]-result_tab$Par_CalledEM_prop[10:12])


mean(result_tab$Par_CalledAR_prop[10:12])



#
mean(result_tab$Par_CalledAR_prop[1:6])
mean(result_tab$Par_CalledAR_prop[7:12])


mean(result_tab$Par_CalledAR_prop[1:3]/result_tab$Par_CalledAR_prop[7:9])

mean(result_tab$SC_CalledEM_prop[7:9]/result_tab$SC_CalledEM_prop[1:3])
mean(result_tab$Par_CalledEM_prop[7:9]/result_tab$Par_CalledEM_prop[1:3])

mean((result_tab$Par_CalledEM_prop[1:6] *result_tab$Total_Par_SNPs[1:6] + result_tab$SC_CalledEM_prop[1:6] *result_tab$total_SC_SNPs[1:6])/(result_tab$total_SC_SNPs[1:6]+result_tab$Total_Par_SNPs[1:6]))

mean((result_tab$Par_CalledEM_prop[7:12] *result_tab$Total_Par_SNPs[7:12] + result_tab$SC_CalledEM_prop[7:12] *result_tab$total_SC_SNPs[7:12])/(result_tab$total_SC_SNPs[7:12]+result_tab$Total_Par_SNPs[7:12]))


mean(((result_tab$Par_CalledEM_prop[7:9] *result_tab$Total_Par_SNPs[7:9] + result_tab$SC_CalledEM_prop[7:9] *result_tab$total_SC_SNPs[7:9])/(result_tab$total_SC_SNPs[7:9]+result_tab$Total_Par_SNPs[7:9]))/((result_tab$Par_CalledEM_prop[1:3] *result_tab$Total_Par_SNPs[1:3] + result_tab$SC_CalledEM_prop[1:3] *result_tab$total_SC_SNPs[1:3])/(result_tab$total_SC_SNPs[1:3]+result_tab$Total_Par_SNPs[1:3])))



mean(result_tab$Par_CalledEM_prop[1:6]+result_tab$Par_CalledAR_prop[1:6])
mean(result_tab$Par_CalledwrongEM_prop[1:6])
mean(1- result_tab$Par_CalledEM_prop[1:6]-result_tab$Par_CalledAR_prop[1:6]-result_tab$Par_CalledwrongEM_prop[1:6])

mean(result_tab$SC_CalledEM_prop[1:6])
mean(result_tab$SC_CalledwrongEM_prop[1:6]+result_tab$SC_CalledwrongAR_prop[1:6])
mean(1- result_tab$SC_CalledEM_prop[1:6] - result_tab$SC_CalledwrongEM_prop[1:6]-result_tab$SC_CalledwrongAR_prop[1:6])

#
mean(result_tab$Par_CalledEM_prop[7:12]+result_tab$Par_CalledAR_prop[7:12])
mean(result_tab$Par_CalledwrongEM_prop[7:12])
mean(1- result_tab$Par_CalledEM_prop[7:12]-result_tab$Par_CalledAR_prop[7:12]-result_tab$Par_CalledwrongEM_prop[7:12])

mean(result_tab$SC_CalledEM_prop[7:12])
mean(result_tab$SC_CalledwrongEM_prop[7:12]+result_tab$SC_CalledwrongAR_prop[7:12])
mean(1- result_tab$SC_CalledEM_prop[7:12] - result_tab$SC_CalledwrongEM_prop[7:12]-result_tab$SC_CalledwrongAR_prop[7:12])


mean(result_tab$Par_CalledEM_prop[7:12]+result_tab$Par_CalledAR_prop[7:12])


mean(result_tab$Par_CalledEM_prop[7:12])

mean(1/((result_tab$Par_CalledEM_prop[7:12]+result_tab$Par_CalledAR_prop[7:12])/(result_tab$Par_CalledEM_prop[1:6]+result_tab$Par_CalledAR_prop[1:6])))
mean(result_tab$SC_CalledEM_prop[7:12]/result_tab$SC_CalledEM_prop[1:6])
mean(1/((1- result_tab$SC_CalledEM_prop[7:12] - result_tab$SC_CalledwrongEM_prop[7:12]-result_tab$SC_CalledwrongAR_prop[7:12])/(1- result_tab$SC_CalledEM_prop[1:6] - result_tab$SC_CalledwrongEM_prop[1:6]-result_tab$SC_CalledwrongAR_prop[1:6])))

###%correct indentified SNPs
mean(result_tab$SC_correctlyCall_prop[1:3])
mean(result_tab$SC_notCalled_prop[1:3])

mean(result_tab$Par_correctlyCall_prop[1:3])
mean(result_tab$Par_notCalled_prop[1:3])
mean(result_tab$Called_prop[1:3])

mean(result_tab$SC_correctlyCall_prop[4:6])
mean(result_tab$SC_notCalled_prop[4:6])

mean(result_tab$Par_correctlyCall_prop[4:6])
mean(result_tab$Par_notCalled_prop[4:6])
mean(result_tab$Called_prop[4:6])

mean(result_tab$SC_correctlyCall_prop[7:9])
mean(result_tab$SC_notCalled_prop[7:9])
mean(result_tab$Called_prop[7:9])

mean(result_tab$Par_correctlyCall_prop[7:9])
mean(result_tab$Par_notCalled_prop[7:9])
mean(result_tab$Called_prop[7:9])

mean(result_tab$SC_correctlyCall_prop[10:12])
mean(result_tab$SC_notCalled_prop[10:12])
  
mean(result_tab$Par_correctlyCall_prop[10:12])
mean(result_tab$Par_notCalled_prop[10:12])
mean(result_tab$Called_prop[10:12])

result_tab$Called_prop[7:9]-result_tab$Called_prop[10:12]
mean(result_tab$Called_prop[7:9]-result_tab$Called_prop[10:12])


mean(1/(result_tab$Called_prop[c(4:6, 10:12)])/mean(result_tab$Called_prop[c(1:3,7:9)]))
mean(1/(result_tab$Called_prop[c(4:6)])/mean(result_tab$Called_prop[c(1:3)]))
mean(1/(result_tab$Called_prop[c(10:12)])/mean(result_tab$Called_prop[c(7:9)]))

mean(1/(result_tab$Called_prop[c(7:12)])/mean(result_tab$Called_prop[c(1:6)]))
mean(1/result_tab$Called_prop[c(7:9)])/mean(result_tab$Called_prop[c(1:3)])
mean(1/result_tab$Called_prop[c(10:12)])/mean(result_tab$Called_prop[c(4:6)])



mean(1/(result_tab$Called_prop[c(7:12)])/mean(result_tab$Called_prop[c(1:3,7:9)]))


mean(1/(result_tab$SC_correctlyCall_prop[c(4:6)])/mean(result_tab$SC_correctlyCall_prop[c(1:3)]))
1-mean(result_tab$SC_correctlyCall_prop[c(10:12)])/mean(result_tab$SC_correctlyCall_prop[c(7:9)])
1-mean(result_tab$Par_correctlyCall_prop[c(4:6)])/mean(result_tab$Par_correctlyCall_prop[c(1:3)])
1-mean(result_tab$Par_correctlyCall_prop[c(10:12)])/mean(result_tab$Par_correctlyCall_prop[c(7:9)])

###

mean(result_tab$Agreement_in_SingleCopyRegions[1:3])
mean(result_tab$Agreement_in_MultiCopyRegions[1:3])

mean(result_tab$Agreement_in_SingleCopyRegions[4:6])
mean(result_tab$Agreement_in_MultiCopyRegions[4:6])



mean(result_tab$Agreement_in_SingleCopyRegions[7:9])
mean(result_tab$Agreement_in_MultiCopyRegions[7:9])

mean(result_tab$Agreement_in_SingleCopyRegions[10:12])
mean(result_tab$Agreement_in_MultiCopyRegions[10:12])




mean(result_tab$ParalogTaggedbyParaMask[1:3]/result_tab$ParalogsWithSNPs[1:3])
mean(result_tab$ParalogTaggedbyParaMask[1:3]/result_tab$ParalogsWithSNPs[1:3])

mean(result_tab$ParalogTaggedbyParaMask[4:6]/result_tab$ParalogsWithSNPs[4:6])
mean(result_tab$ParalogTaggedbyParaMask[4:6]/result_tab$ParalogsWithSNPs[4:6])


mean(result_tab$ParalogTaggedbyParaMask[1:3])
mean(result_tab$ParalogTaggedbyParaMask[4:6]/result_tab$ParalogsWithSNPs[4:6])


mean(result_tab$Agreement_in_MultiCopyRegions[1:3])

####check overlap?


##for plotting later

# SVtab_infer<-SVtab_infer2

SVtab$data <- "real"

SVtab_infer$data  <- "estimated"

SVtab$data<- factor(SVtab$data, levels = c("real", "estimated"))
SVtab_infer$data<- factor(SVtab_infer$data, levels = c("real", "estimated"))


tmp<-overlap$overlap[1]
start<- overlap$pos[1]
extend<-TRUE
overlap_tab<-c()
sum(overlap$overlap==0)
for(i in 2:nrow(overlap)){
  tmp_new<-overlap$overlap[i]
  if(!extend){
    tmp<-overlap$overlap[i]
    start<- overlap$pos[i]
    extend <- TRUE
  }else if(tmp!=tmp_new){
    end <- overlap$pos[i]
    extend<- FALSE
    overlap_tab<- rbind(overlap_tab, c(start,end, tmp))
  }
}
end <- overlap$pos[length(overlap$pos)]
overlap_tab<- rbind(overlap_tab, c(start,end, tmp))

head(overlap_tab)
overlap_tab
overlap_tab <- as.data.frame(overlap_tab)
colnames(overlap_tab) <- c("start", "end", "overlap")
head(overlap_tab)
overlap_tab$data <- "overlap"
overlap_tab$overlap<- as.factor(overlap_tab$overlap)


tmp<-overlap_tab[overlap_tab$overlap==2,]
tmp2<- SVtab[SVtab$type=="par",]
sum(tmp$end-tmp$start)/sum(tmp2$end- tmp2$start)

tmp<-overlap_tab[overlap_tab$overlap==1,]
tmp2<- SVtab[SVtab$type=="sc",]
sum(tmp$end-tmp$start)/sum(tmp2$end- tmp2$start)

overlap_tab$overlap2<-"joint"
overlap_tab$overlap2[overlap_tab$overlap==0 | overlap_tab$overlap==3]<-"disjoint"

SVtab$ypos<-2
overlap_tab$ypos<-1
SVtab_infer$ypos<-0




g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 




SVtab_infer$data <-factor(as.character(SVtab_infer$data), levels=c( "position","real", "overlap", "estimated"))
SVtab$data <-factor(as.character(SVtab$data), levels=c( "position","real", "overlap", "estimated"))
overlap_tab$data <-factor(as.character(overlap_tab$data), levels=c( "position","real", "overlap", "estimated"))
pos1<- data.frame(start=0, end=500000, data="position")
pos1$data <-factor(as.character(pos1$data), levels=c( "position","real", "overlap", "estimated"))
pos2<- data.frame(start=500000, end=1000000, data="position")
pos2$data <-factor(as.character(pos2$data), levels=c( "position","real", "overlap", "estimated"))



SVtab_infer$size<-15
SVtab$size<-15
SVtab_infer$size[SVtab_infer$type==1]<-20
SVtab$size[SVtab$type=="par"]<-20


SVtab_infer_p1<- SVtab_infer[SVtab_infer$Start<=500000,]
SVtab_infer_p1$End[nrow(SVtab_infer_p1)]<-500000
SVtab_infer_p2<- SVtab_infer[SVtab_infer$End>=500000,]
SVtab_infer_p2$Start[1]<-500000
SVtab_infer_p1$type[SVtab_infer_p1$type.0.single.copy.1.multi.copy==0] <- "sc"
SVtab_infer_p1$type[SVtab_infer_p1$type.0.single.copy.1.multi.copy==1] <- "par"
SVtab_infer_p2$type[SVtab_infer_p2$type.0.single.copy.1.multi.copy==0] <- "sc"
SVtab_infer_p2$type[SVtab_infer_p2$type.0.single.copy.1.multi.copy==1] <- "par"


SVtab_real_p1<- SVtab[SVtab$start<=500000,]
SVtab_real_p1$end[nrow(SVtab_real_p1)]<-500000
SVtab_real_p2<- SVtab[SVtab$end>=500000,]
SVtab_real_p2$start[1]<-500000

overlap_tab_p1<- overlap_tab[overlap_tab$start<=500000,]
overlap_tab_p1$end[nrow(overlap_tab_p1)]<-500000
overlap_tab_p2<- overlap_tab[overlap_tab$end>=500000,]
overlap_tab_p2$start[1]<-500000

overlap_tab
SVtab_real_p2$type
SVtab_infer_p2$type


segplot <- ggplot()+
  geom_segment(data = SVtab_infer_p1, mapping = aes(x = Start, y = data, xend=End, yend=data, color=type, size=size))+
  geom_segment(data = SVtab_real_p1, mapping = aes(x = start, y = data, xend=end, yend=data, color=type, size=size))+
  geom_segment(data = overlap_tab_p1[overlap_tab_p1$overlap2=="joint",], mapping = aes(x = start, y = data, xend=end, yend=data, color=overlap2), size=10)+
  geom_segment(data = pos1, aes(x = start, y =data, xend=end, yend=data), color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=500000, xend=500000, color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=0, xend=0, color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=250000, xend=250000, color="black", size=1)+
  annotate(geom = "text", y=0.4, x=500000, color="black", label="500 kbp",size=5)+
  annotate(geom = "text", y=0.4, x=0, color="black", label="0 kbp", size=5)+
  annotate(geom = "text", y=0.4, x=250000, color="black",label="250 kbp", size=5)+
  scale_y_discrete( expand = c(0,0), limits=c("position","estimated", "overlap", "real"), labels = c( "","inferred", "Agreement", "simulated"))+
  scale_x_continuous(expand = c(0,0), limits = c(-13000,513000), breaks=seq(0,1000000, by=250000), labels=seq(0,1000, by=250))+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.text = element_text(size=16), legend.position = "top", aspect.ratio = 1, plot.margin=unit(c(-20,0,-20,0), "cm"), panel.margin=unit(0,"null"), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.key.width = unit(5, 'cm'), legend.key.height  = unit(5, 'cm'), legend.key = element_rect(fill = "white"), legend.text = element_text(size=20))+
  labs(x="", y="", color="")+
  scale_color_manual(limits=c("sc", "par"),breaks = c("sc", "par", "joint", "disjoint"),labels = c("single copy", "paralogs", "joint", "disjoint"),values=c("blue3","brown2","honeydew4", "black"))+
  guides(color=guide_legend(title = ""), size="none")

segplot_legend<-g_legend(segplot)
plot_grid(segplot_legend)

 overlap_plot<- ggplot()+
  geom_segment(data = SVtab_infer_p1[SVtab_infer_p1$type=="sc",], mapping = aes(x = Start, y = data, xend=End, yend=data, color=type), size=4,position = position_nudge(y=-5))+
  geom_segment(data = SVtab_infer_p1[SVtab_infer_p1$type=="par",], mapping = aes(x = Start, y = data, xend=End, yend=data, color=type), size=25, position = position_nudge(y=-5))+
  geom_segment(data = SVtab_real_p1[SVtab_real_p1$type=="sc",], mapping = aes(x = start, y = data, xend=end, yend=data, color=type),size=4,position = position_nudge(y=5))+
  geom_segment(data = SVtab_real_p1[SVtab_real_p1$type=="par",], mapping = aes(x = start, y = data, xend=end, yend=data, color=type),size=25,position = position_nudge(y=5))+
  # geom_segment(data = SVtab_real_p1, mapping = aes(x = start, y = data, xend=end, yend=data, color=type, size=size),size=3)+
  geom_segment(data = overlap_tab_p1[overlap_tab_p1$overlap2=="joint",], mapping = aes(x = start, y = data, xend=end, yend=data, color=overlap2), size=25, )+
  geom_segment(data = pos1, aes(x = start, y =data, xend=end, yend=data), color="black", size=1, position = position_nudge(y=-10))+
  annotate(geom = "segment", y=-9, yend=-10, x=500000, xend=500000, color="black", size=1)+
  annotate(geom = "segment", y=-9, yend=-10, x=0, xend=0, color="black", size=1)+
  annotate(geom = "segment", y=-9, yend=-10, x=250000, xend=250000, color="black", size=1)+
  annotate(geom = "text", y=-12, x=500000, color="black", label="500 kbp",size=13)+
  annotate(geom = "text", y=-12, x=0, color="black", label="0 kbp", size=13)+
  annotate(geom = "text", y=-12, x=250000, color="black",label="250 kbp", size=13)+
  annotate(geom = "text", y=9.25, x=-50000, color="black", label="Simulated",size=13)+
  annotate(geom = "text", y=3, x=-50000, color="black", label="Agreement", size=13)+
  annotate(geom = "text", y=-2.75, x=-50000, color="black",label="Inferred", size=13)+
  scale_y_discrete( expand = c(6,0), limits=c("position","estimated", "overlap", "real"), labels = c( "","inferred", "agreement", "simulated"))+
  scale_x_continuous(expand = c(0,0), limits = c(-100000,550000), breaks=seq(0,1000000, by=250000), labels=seq(0,1000, by=250))+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.text = element_text(size=16), legend.position = "", panel.margin=unit(0,"null"),plot.margin = margin(t=20, r=1, l=1, b=1), axis.ticks.y = element_blank(), axis.ticks = element_line(size=5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank())+
  labs(x="", y="", color="")+
  scale_color_manual(limits=c("sc", "par"),breaks = c("sc", "par", "joint", "disjoint"),labels = c("single copy", "paralogs", "joint", "disjoint"),values=c("blue3","brown2","honeydew4", "black"))+
  guides(color=guide_legend(title = ""), size="none")




ggplot()+
  geom_segment(data = SVtab_infer_p2[SVtab_infer_p2$type=="sc",], mapping = aes(x = Start, y = data, xend=End, yend=data, color=type), size=2)+
  geom_segment(data = SVtab_infer_p2[SVtab_infer_p2$type=="par",], mapping = aes(x = Start, y = data, xend=End, yend=data, color=type), size=10)+
  geom_segment(data = SVtab_real_p2[SVtab_real_p2$type=="sc",], mapping = aes(x = start, y = data, xend=end, yend=data, color=type),size=2)+
  geom_segment(data = SVtab_real_p2[SVtab_real_p2$type=="par",], mapping = aes(x = start, y = data, xend=end, yend=data, color=type),size=10)+
  geom_segment(data = overlap_tab_p2[overlap_tab_p2$overlap2=="joint",], mapping = aes(x = start, y = data, xend=end, yend=data, color=overlap2), size=8)+
  geom_segment(data = pos2, aes(x = start, y =data, xend=end, yend=data), color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=1000000, xend=1000000, color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=500000, xend=500000, color="black", size=1)+
  annotate(geom = "segment", y=1, yend=0.7, x=750000, xend=750000, color="black", size=1)+
  annotate(geom = "text", y=0.4, x=1000000, color="black", label="1000 kbp",size=5)+
  annotate(geom = "text", y=0.4, x=500000, color="black", label="500 kbp",size=5)+
  annotate(geom = "text", y=0.4, x=750000, color="black",label="750 kbp",size=5)+
  scale_y_discrete( expand = c(6,0), limits=c("position","estimated", "overlap", "real"), labels = c( "","inferred", "agreement", "simulated"))+
  scale_x_continuous(expand = c(0,0), limits = c(487000,1013000), breaks=seq(0,1000000, by=250000), labels=seq(0,1000, by=250))+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.text = element_text(size=16), legend.position = "top", aspect.ratio = 1, plot.margin=unit(c(-20,0,-20,0), "cm"), panel.margin=unit(0,"null"), axis.ticks.y = element_blank(),  axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  labs(x="", y="", color="")+
  scale_color_manual(breaks = c("sc", "par", "joint", "disjoint"),labels = c("single copy", "paralogs", "joint", "disjoint"),values=c("blue3","brown2","honeydew4", "black"))+
  guides(size="none", color="none")























# ###check single classes 
# clusters<- read.table(file = "/home/btjeng/Data/SeDuS/Sim_final2/HW_10PercentParalogs/ParaMask_runs/Sim0.1HW_rep1EMresults.clusters.txt", header = T, sep = "\t")
# 
# head(clusters)
# 
# ###why????????????????????????????????????????????
# length(unique(clusters$cluster))
# sum(SVtab_infer$type.0.single.copy.1.multi.copy==1)
# tail(SVtab_infer)
# 
# clusters$cluster<- factor(clusters$cluster, levels = 1:318)
# sum(tabulate(clusters$cluster)==1)/length(levels(clusters$cluster))
# 
# SVtab_infer$nSeeds <-0
# for(i in 1:nrow(SVtab_infer)){
#   SVtab_infer$nSeeds[i] <- sum(clusters$position>=SVtab_infer$Start[i] & clusters$position<=SVtab_infer$End[i])
# }
# sum(SVtab_infer$nSeeds==1)
# 
# 
# overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000),inferred_cluster=rep(0, 1000000), nSeeds=rep(0, 1000000))
# head(overlap)
# c<-1
# for(i in 1:nrow(SVtab_infer)){
#   if(SVtab_infer$type[i]==1){
#     overlap$inferred[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- 1
#     overlap$nSeeds[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- SVtab_infer$nSeeds[i]
#     overlap$inferred_cluster[overlap$pos >= SVtab_infer$Start[i] & overlap$pos <=SVtab_infer$End[i]] <- c
#     c <- c + 1
#   }
# }
# 
# 
# for(i in 1:nrow(SVtab)){
#   if(SVtab$type[i]=="par"){
#     overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <- 1
#   }
# }
# head(overlap)
# tmp1<-overlap[overlap$real==0,]
# tmp2<-overlap[overlap$real==1,]
# 
# unique(tmp1$inferred_cluster)
# unique(tmp2$inferred_cluster)
# 
# unique(tmp1$nSeeds)
# unique(tmp2$nSeeds)
# tabulate(factor(tmp2$nSeeds, levels = unique(tmp2$nSeeds)))
# 
# tmp1.1<-tmp1[tmp1$nSeeds==1,]
# tmp2.1<-tmp2[tmp2$nSeeds==1,]
# 
# length(unique(tmp1.1$inferred_cluster))
# length(unique(tmp2.1$inferred_cluster))
# 
# unique(overlap$inferred_cluster)
# unique(tmp1.1$inferred_cluster)
# head(SVtab_infer)
# 
# SVtab_infer$cluster<-0
# SVtab_infer$cluster[SVtab_infer$type.0.single.copy.1.multi.copy==1] <- 1:sum(SVtab_infer$type.0.single.copy.1.multi.copy==1)
# tail(SVtab_infer)
# sum(SVtab_infer$nSeeds==1)
# SVtab_infer2<-SVtab_infer
# SVtab_infer2$type[SVtab_infer2$nSeeds==1]<-0
# SVtab_infer2$type.0.single.copy.1.multi.copy[SVtab_infer2$nSeeds==1]<-0
# 
# head(clusters)
# Sim$cluster <- 0
# 
# i<-1
# for(i in 1:nrow(clusters)){
#   Sim$cluster[Sim$Position==clusters$position[i]] <-clusters$cluster[i]
# }
# 
# sum(Sim$cluster==13)
# Sim$finalClass[Sim$cluster==13]
# Sim$finalClass[Sim$cluster %in% unique(SVtab_infer$cluster[SVtab_infer$nSeeds==1])]<-0
# 
# 
# ###testing overlap with filtered
# overlap<- data.frame(pos=c(1:1000000), real=rep(0, 1000000),inferred=rep(0, 1000000))
# head(overlap)
# for(i in 1:nrow(SVtab_infer2)){
#   if(SVtab_infer2$type[i]==1){
#     overlap$inferred[overlap$pos >= SVtab_infer2$Start[i] & overlap$pos <=SVtab_infer2$End[i]] <- 1
#   }
# }
# 
# for(i in 1:nrow(SVtab)){
#   if(SVtab_real$type[i]=="par"){
#     overlap$real[overlap$pos >=SVtab$start[i] & overlap$pos <=SVtab$end[i]] <- 1
#   }
# }
# 
# head(overlap)
# overlap$overlap<-0
# overlap$overlap[overlap$real==0 & overlap$inferred==0]<- 1
# overlap$overlap[overlap$real==1 & overlap$inferred==1]<- 2
# overlap$overlap[overlap$real==1 & overlap$inferred==0]<- 3
# 
# sc_agreement <- sum(overlap$overlap==1)/sum(overlap$real==0)
# sc_discrepency <- sum(overlap$overlap==0)/sum(overlap$real==0)
# 
# par_agreement <- sum(overlap$overlap==2)/sum(overlap$real==1)
# par_discrepency <- sum(overlap$overlap==3)/sum(overlap$real==1)
# 
# overall_agreement <- sum(overlap$overlap==1 | overlap$overlap==2)/nrow(overlap)


## make plot instead of table

ptab_HW <- rbind(cbind(result_tab$Paralog_prop[result_tab$Fis==0], result_tab$Par_correctlyCall_prop[result_tab$Fis==0], "par"),cbind(result_tab$Paralog_prop[result_tab$Fis==0], result_tab$SC_correctlyCall_prop[result_tab$Fis==0], "sc"))

ptab_HW <- as.data.frame(ptab_HW)
colnames(ptab_HW) <- c("Prop", "Correct", "type")
ptab_HW$Prop <- as.numeric(ptab_HW$Prop)
ptab_HW$Correct <- as.numeric(ptab_HW$Correct)
ptab_HW$Replicate <- rep(1:3, 4)
ptab_HW$Replicate<- as.factor(ptab_HW$Replicate)

ptab_F <- rbind(cbind(result_tab$Paralog_prop[result_tab$Fis==0.9], result_tab$Par_correctlyCall_prop[result_tab$Fis==0.9], "par"),cbind(result_tab$Paralog_prop[result_tab$Fis==0.9], result_tab$SC_correctlyCall_prop[result_tab$Fis==0.9], "sc"))

ptab_F <- as.data.frame(ptab_F)
colnames(ptab_F) <- c("Prop", "Correct", "type")
ptab_F$Prop <- as.numeric(ptab_F$Prop)
ptab_F$Correct <- as.numeric(ptab_F$Correct)
ptab_F$Replicate <- rep(1:3, 4)
ptab_F$Replicate<- as.factor(ptab_F$Replicate)


plot_HW<-ggplot(data = ptab_HW, aes(x = Prop, y=Correct, color=type, shape=Replicate))+
  geom_point(size=6,stroke=2)+
  scale_x_continuous(breaks = c(0.1,0.5), labels=c(10,50), limits = c(0.05,0.55), expand=c(0,0))+
  labs(y= "Recall (%)", x="Proportion of duplications (%)")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x =  element_text(size=40, vjust = -0.55), axis.title.y = element_text(size=40, vjust = 1), axis.text = element_text(size=40), legend.title = element_blank(), legend.text = element_text(size=14), legend.position = "", plot.margin = margin(t=90, r=60, l=1, b=60))+
  scale_y_continuous(breaks = seq(0.88, 1, by=0.06), limits = c(0.85,1.01), expand=c(0,0), labels = seq(88,100, by=6))+
  geom_segment(x=0.05, xend=0.05, y=0.88, yend=1, color="black")+
  geom_segment(x=0.1, xend=0.5, y=0.85, yend=0.85, color="black")+
  # scale_color_brewer(palette = "Dark2",breaks=c("Collapsed", "ParaOut", "SC"), labels=c("collapsed", "masked", "simulated single copy"))+
  scale_color_manual(values=c("brown2", "blue3"),breaks=c("par", "sc"), labels=c("multicopy", "single-copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))

plot_FI<-ggplot(data = ptab_F, aes(x = Prop, y=Correct, color=type, shape=Replicate))+
  geom_point(size=6,stroke=2)+
  scale_x_continuous(breaks = c(0.1,0.5), labels=c(10,50), limits = c(0.05,0.55), expand=c(0,0))+
  labs(y= "Proportion of correctly identified SNPs (%)", x="Proportion of duplications (%)")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x =  element_text(size=40, vjust = -0.55), axis.title.y = element_blank(), axis.text.x = element_text(size=40), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size=14), legend.position = "", plot.margin = margin(t=90, r=60, l=1, b=60))+
  scale_y_continuous(breaks = seq(0.88, 1, by=0.06), limits = c(0.85,1.01), expand=c(0,0), labels = seq(88,100, by=6))+
  geom_segment(x=0.1, xend=0.5, y=0.85, yend=0.85, color="black")+
  scale_color_manual(values=c("brown2", "blue3"),breaks=c("par", "sc"), labels=c("multicopy", "single-copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))

l_plot<-ggplot(data = ptab_F, aes(x = Prop, y=Correct, color=type, shape=Replicate))+
  geom_point(size=6,stroke=2)+
  scale_x_continuous(breaks = c(0.1,0.5), labels=c(10,50), limits = c(0.05,0.55), expand=c(0,0))+
  labs(y= "Proportion of correctly identified SNPs (%)", x="Proportion of duplications (%)")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x =  element_text(size=40, vjust = -0.55), axis.title.y = element_blank(), axis.text.x = element_text(size=40), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size=40), plot.margin = margin(t=90, r=60, l=1, b=60),legend.spacing.y = unit(1.0,"cm"), legend.key.size = unit(1.0,"cm"), legend.key = element_rect(fill=NA, color = "white"))+
  scale_y_continuous(breaks = seq(0.88, 1, by=0.06), limits = c(0.85,1.01), expand=c(0,0), labels = seq(88,100, by=6))+
  geom_segment(x=0.1, xend=0.5, y=0.85, yend=0.85, color="black")+
  scale_color_manual(values=c("brown2", "blue3"),breaks=c("par", "sc"), labels=c("multicopy", "single-copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))+
  guides(shape=F, color=guide_legend(byrow = TRUE,override.aes = list(size = 10)))

plot_grid(plotlist=list(plot_HW, plot_FI), ncol=2, align='hv')

ptab_HW <- rbind(cbind(result_tab$Paralog_prop[result_tab$Fis==0], result_tab$Par_correctlyCall_prop[result_tab$Fis==0], "par"),cbind(result_tab$Paralog_prop[result_tab$Fis==0], result_tab$SC_correctlyCall_prop[result_tab$Fis==0], "sc"))


g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 


legend<-g_legend(l_plot)
plot_grid(legend)
recall_plot<-plot_grid(plotlist=list(plot_HW, plot_FI), ncol=2, align='hv')
plot_grid(plotlist=list(recall_plot, legend), ncol=2, rel_widths = c(0.83,0.17))



##seeds contributions

HW_em_par_0.1 <- result_tab$Par_CalledEM_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0]
HW_em_sc_0.1 <- result_tab$SC_CalledEM_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0]
HW_em_par_0.5 <- result_tab$Par_CalledEM_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0]
HW_em_sc_0.5 <- result_tab$SC_CalledEM_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0]

HW_AR_par_0.1 <- result_tab$Par_CalledAR_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0]
HW_AR_par_0.5 <- result_tab$Par_CalledAR_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0]

HW_total_par_0.1 <- result_tab$Par_correctlyCall_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0]
HW_total_sc_0.1 <- result_tab$SC_correctlyCall_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0]
HW_total_par_0.5 <- result_tab$Par_correctlyCall_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0]
HW_total_sc_0.5 <- result_tab$SC_correctlyCall_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0]

#
F_em_par_0.1 <- result_tab$Par_CalledEM_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0.9]
F_em_sc_0.1 <- result_tab$SC_CalledEM_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0.9]
F_em_par_0.5 <- result_tab$Par_CalledEM_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0.9]
F_em_sc_0.5 <- result_tab$SC_CalledEM_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0.9]

F_AR_par_0.1 <- result_tab$Par_CalledAR_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0.9]
F_AR_par_0.5 <- result_tab$Par_CalledAR_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0.9]

F_total_par_0.1 <- result_tab$Par_correctlyCall_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0.9]
F_total_sc_0.1 <- result_tab$SC_correctlyCall_prop[result_tab$Paralog_prop==0.1 & result_tab$Fis==0.9]
F_total_par_0.5 <- result_tab$Par_correctlyCall_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0.9]
F_total_sc_0.5 <- result_tab$SC_correctlyCall_prop[result_tab$Paralog_prop==0.5 & result_tab$Fis==0.9]


stab <- rbind(cbind(0.1,0, mean(HW_total_par_0.1), sd(HW_total_par_0.1), "par", "Clustering"),
                 cbind(0.1,0, mean(HW_em_par_0.1), sd(HW_em_par_0.1), "par", "EM"),
                 cbind(0.1,0, mean(HW_AR_par_0.1), sd(HW_AR_par_0.1), "par", "AR"),
                 cbind(0.1,0, mean(HW_total_sc_0.1), sd(HW_total_sc_0.1), "sc", "Clustering"),
                 cbind(0.1,0, mean(HW_em_sc_0.1), sd(HW_em_sc_0.1), "sc", "EM"),
                 
                 cbind(0.5,0, mean(HW_total_par_0.5), sd(HW_total_par_0.5), "par", "Clustering"),
                 cbind(0.5,0, mean(HW_em_par_0.5), sd(HW_em_par_0.5), "par", "EM"),
                 cbind(0.5,0, mean(HW_AR_par_0.5), sd(HW_AR_par_0.5), "par", "AR"),
                 cbind(0.5,0, mean(HW_total_sc_0.5), sd(HW_total_sc_0.5), "sc", "Clustering"),
                 cbind(0.5,0, mean(HW_em_sc_0.5), sd(HW_em_sc_0.5), "sc", "EM"),
                 
                 cbind(0.1,0.9, mean(F_total_par_0.1), sd(F_total_par_0.1), "par", "Clustering"),
                 cbind(0.1,0.9, mean(F_em_par_0.1), sd(F_em_par_0.1), "par", "EM"),
                 cbind(0.1,0.9, mean(F_AR_par_0.1), sd(F_AR_par_0.1), "par", "AR"),
                 cbind(0.1,0.9, mean(F_total_sc_0.1), sd(F_total_sc_0.1), "sc", "Clustering"),
                 cbind(0.1,0.9, mean(F_em_sc_0.1), sd(F_em_sc_0.1), "sc", "EM"),
                 
                 cbind(0.5,0.9, mean(F_total_par_0.5), sd(F_total_par_0.5), "par", "Clustering"),
                 cbind(0.5,0.9, mean(F_em_par_0.5), sd(F_em_par_0.5), "par", "EM"),
                 cbind(0.5,0.9, mean(F_AR_par_0.5), sd(F_AR_par_0.5), "par", "AR"),
                 cbind(0.5,0.9, mean(F_total_sc_0.5), sd(F_total_sc_0.5), "sc", "Clustering"),
                 cbind(0.5,0.9, mean(F_em_sc_0.5), sd(F_em_sc_0.5), "sc", "EM")
                 
)
stab <- as.data.frame(stab)
colnames(stab) <- c("Prop","Fis", "Recall", "sd", "type", "classifier")

stab$Prop <- as.numeric(stab$Prop)
stab$Recall <- as.numeric(stab$Recall)

stab$plotfac <- paste(stab$Prop,stab$Fis, sep = "_")
stab$plotfac <- factor(stab$plotfac, levels= unique(stab$plotfac)[c(4:1)])

stab$classifier <- factor(stab$classifier, levels = c("Clustering", "AR", "EM"))
stab$group <- paste(stab$type,stab$classifier, sep = "_")
stab$group <- factor(stab$group, levels = unique(stab$group)[length(unique(stab$group)):1])
stab$upper <-stab$Recall+as.numeric(stab$sd)
stab$lower <-stab$Recall-as.numeric(stab$sd)





# write.table(x = stab, file = "~/Data/Paramask/Sim_final/stab.txt", sep="\t", col.names = T)
stab<- read.table("~/Data/Paramask/Sim_final/stab.txt", header = T, sep = "\t")

stab$plotfac


guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- y$label_trans(trained$key$.label)
  trained
}
ggplot(data = stab, aes(x = plotfac, y=Recall, group=group))+
  labs(x="Proportion of duplications (%)", y="Recall (%)")+
  geom_bar(aes(fill=group , alpha=group),stat="identity", position = "dodge", width = 0.5,color="black", size=1.1)+
  geom_errorbar(aes(ymin=lower, ymax=upper,group=group), width=0.2, position=position_dodge(width=0.5), size=1.1)+
  scale_x_discrete(limits=c(levels(stab$plotfac)), breaks=levels(stab$plotfac), labels= c(50,10,50,10))+
  scale_y_continuous(position = "right",breaks=seq(0,1, by=0.5), labels = seq(0,100, by=50), limits=c(-0.025,1.25), expand=c(0,0))+
  scale_fill_manual(values=c("brown2", "brown2", "brown2", "blue3","blue3"),breaks=c("par_Clustering","par_EM", "par_AR","sc_Clustering","sc_EM" ), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_alpha_manual(breaks = c("par_Clustering","par_EM", "par_AR","sc_Clustering","sc_EM" ), values = c(1, 0.6, 0.2, 1, 0.6),labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  geom_segment(data = stab,aes(x="0.1_0", xend="0.1_0", y=0, yend=1), color="black", position = position_nudge(x = 0.6))+
  geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=-0.025, yend=-0.025), color="black")+
  geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=-0.025, yend=-0.025), color="black")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x.top  =  element_text(size=40, vjust = 4, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 9) , axis.text = element_text(size=40), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "bottom", plot.margin = margin(t=30, r=30, l=30, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))+
  coord_flip()+
  guides(
    fill = guide_legend(override.aes = list( linetype = 0), nrow=3, byrow = T)
  )

splot<- ggplot(data = stab, aes(x = plotfac, y=Recall, group=group))+
  labs(x="Proportion of duplications (%)", y="Recall")+
  geom_bar(aes(fill=group , alpha=group),stat="identity", position = "dodge", width = 0.5,color="black", size=1.1)+
  geom_errorbar(aes(ymin=lower, ymax=upper,group=group), width=0.2, position=position_dodge(width=0.5), size=1.1)+
  scale_x_discrete(limits=c(levels(stab$plotfac)), breaks=levels(stab$plotfac), labels= c(50,10,50,10))+
  scale_y_continuous(position = "right",breaks=seq(0,1, by=0.5), limits=c(-0.025,1.025), expand=c(0,0))+
  scale_fill_manual(values=c("brown2", "brown2", "brown2", "blue3","blue3"),breaks=c("par_Clustering","par_EM", "par_AR","sc_Clustering","sc_EM" ), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_alpha_manual(breaks = c("par_Clustering","par_EM", "par_AR","sc_Clustering","sc_EM" ), values = c(1, 0.6, 0.2, 1, 0.6),labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  geom_segment(data = stab,aes(x="0.1_0", xend="0.1_0", y=0, yend=1), color="black", position = position_nudge(x = 0.6))+
  geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=-0.025, yend=-0.025), color="black")+
  geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=-0.025, yend=-0.025), color="black")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x.top  =  element_text(size=40, vjust = 4, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 10) , axis.text = element_text(size=40), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))+
  coord_flip()+
  guides(
    fill = guide_legend(override.aes = list( linetype = 0), nrow=3, byrow = T)
  )
###


# plot_grid(plotlist=list(splot_HW, splot_F), ncol=2)


###EM plot 
em_plot1 <- plot_grid(plotlist=list(em_plot_HW, em_plot_F), ncol=1, align = "vh" )

y.grob <- textGrob("Heterozygote frequency", 
                   gp=gpar( fontsize=40), rot=90, vjust = 1)

#add to plot

em_plot<- grid.arrange(arrangeGrob(em_plot1, left = y.grob))


em_s_plot <- plot_grid(plotlist=list(em_plot,splot), ncol=2, rel_widths = c(1,1))
em_s_plot
overlap_plot
all_plot <- plot_grid(plotlist=list(em_s_plot, overlap_plot), nrow =2, rel_heights = c(2,1))
all_plot
