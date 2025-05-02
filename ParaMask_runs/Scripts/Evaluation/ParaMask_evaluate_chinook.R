het_chinook<-read.table(file="ParaMask_Chinook3_EMresults_CNstatus.het", header = T, sep = "\t")
head(het_chinook)

het<-read.table(file="ParaMask_Chinook3_EMresults.het", header = T, sep = "\t")

sum(het$EM_class==1 &het$allele.deviation.seed==0)/nrow(het)
sum(het$EM_class==2 |het$allele.deviation.seed==1)/nrow(het)
sum(het$EM_class==0 &het$allele.deviation.seed==0)/nrow(het)


#without uncertain
het_chinook_noUncertain <- het_chinook[!(het_chinook$EM_class==1 &het_chinook$allele.deviation.seed==0), ]
nrow(het_chinook_noUncertain)
het_chinook_noUncertain$ParaMask_final<-1
het_chinook_noUncertain$ParaMask_final[het_chinook_noUncertain$allele.deviation.seed==1 | het_chinook_noUncertain$EM_class==2] <- 2


sum(het_chinook_noUncertain$ParaMask_final==het_chinook_noUncertain$CN_status)/nrow(het_chinook_noUncertain)
sum(het_chinook_noUncertain$ParaMask_final==het_chinook_noUncertain$CN_status)/nrow(het_chinook_noUncertain)

het_chinook_noUncertain_par <- het_chinook_noUncertain[het_chinook_noUncertain$CN_status==2,]
sum(het_chinook_noUncertain_par$ParaMask_final==het_chinook_noUncertain_par$CN_status)/nrow(het_chinook_noUncertain_par)
sum(het_chinook_noUncertain_par$ParaMask_final==1)/nrow(het_chinook_noUncertain_par)


het_chinook_noUncertain_sc <- het_chinook_noUncertain[het_chinook_noUncertain$CN_status==1,]
sum(het_chinook_noUncertain_sc$ParaMask_final==het_chinook_noUncertain_sc$CN_status)/nrow(het_chinook_noUncertain_sc)
sum(het_chinook_noUncertain_sc$ParaMask_final==2)/nrow(het_chinook_noUncertain_sc)


sum(het_chinook_noUncertain$ParaMask_final==2 &het_chinook_noUncertain$CN==2)/sum(het_chinook_noUncertain$ParaMask_final==2)
sum(het_chinook_noUncertain$ParaMask_final==2 &het_chinook_noUncertain$CN==1)/sum(het_chinook_noUncertain$ParaMask_final==2)

sum(het_chinook_noUncertain$ParaMask_final==1 &het_chinook_noUncertain$CN==1)/sum(het_chinook_noUncertain$ParaMask_final==1)


#

het_chinook$ParaMask_final<-1
het_chinook$ParaMask_final[het_chinook$allele.deviation.seed==1 | het_chinook$EM_class==2] <- 2
het_chinook$ParaMask_final[ het_chinook$EM_class==1] <- 3

sum(het_chinook$CN_status==het_chinook$ParaMask_final)/nrow(het_chinook)
het_chinook_par <- het_chinook[het_chinook$CN_status==2,]
sum(het_chinook_par$ParaMask_final==het_chinook_par$CN_status)/nrow(het_chinook_par)
sum(het_chinook_par$ParaMask_final==3)/nrow(het_chinook_par)
sum(het_chinook_par$ParaMask_final==1)/nrow(het_chinook_par)

het_chinook_sc <- het_chinook[het_chinook$CN_status==1,]
sum(het_chinook_sc$ParaMask_final==het_chinook_sc$CN_status)/nrow(het_chinook_sc)
sum(het_chinook_sc$ParaMask_final==3)/nrow(het_chinook_sc)
sum(het_chinook_sc$ParaMask_final==2)/nrow(het_chinook_sc)

sum(het_chinook_par$ParaMask_final!=het_chinook_par$CN_status)/nrow(het_chinook_par)


###

nrow(het_chinook)
het_chinook$ParaMask_final<-1
het_chinook$ParaMask_final[het_chinook$allele.deviation.seed==1 | het_chinook$EM_class==2] <- 2


sum(het_chinook$ParaMask_final==het_chinook$CN_status)/nrow(het_chinook)

het_chinook_par <- het_chinook[het_chinook$CN_status==2,]
sum(het_chinook_par$ParaMask_final==het_chinook_par$CN_status)/nrow(het_chinook_par)
sum(het_chinook_par$ParaMask_final!=het_chinook_par$CN_status)/nrow(het_chinook_par)

het_chinook_sc <- het_chinook[het_chinook$CN_status==1,]
sum(het_chinook_sc$ParaMask_final==het_chinook_sc$CN_status)/nrow(het_chinook_sc)
sum(het_chinook_sc$ParaMask_final!=het_chinook_sc$CN_status)/nrow(het_chinook_sc)


sum(het_chinook$ParaMask_final==2 &het_chinook$CN==2)/sum(het_chinook$ParaMask_final==2)
sum(het_chinook$ParaMask_final==2 &het_chinook$CN==1)/sum(het_chinook$ParaMask_final==2)

sum(het_chinook$ParaMask_final==1 &het_chinook$CN==1)/sum(het_chinook$ParaMask_final==1)


head(het_chinook)
het_chinook$TPN <- 0
het_chinook$TPN[het_chinook$ParaMask_final==1 &het_chinook$CN==1] <-0
het_chinook$TPN[het_chinook$ParaMask_final==1 &het_chinook$CN==2]<-1
het_chinook$TPN[het_chinook$ParaMask_final==2 &het_chinook$CN==1] <-2
het_chinook$TPN[het_chinook$ParaMask_final==2 &het_chinook$CN==2] <-3

het_chinook <- as.data.frame(het_chinook)
het_chinook$EM_class<- factor(het_chinook$EM_class)
het_chinook$CN_status<- factor(het_chinook$CN_status)
het_chinook$TPN<- factor(het_chinook$TPN)

library(ggplot2)
ggplot(het_chinook)+
  geom_point(data = het_chinook[het_chinook$CN_status==1,],aes(x=Heterozygous.geno.freq, y=Het.allele.ratio, color=TPN, shape=TPN, fill=TPN),stroke=1.2, alpha=0.2,size=3)+
  geom_point(data = het_chinook[het_chinook$CN_status==2, ],aes(x=Heterozygous.geno.freq, y=Het.allele.ratio, color=TPN, shape=TPN, fill=TPN),stroke=1.2, alpha=0.2,size=3)+
  scale_x_continuous(limits=c(-0.05,1.05), breaks = seq(0,1,by =0.25), expand=c(0,0))+
  labs(x="Heterozygote frequency", y="Read ratio deviation (D)")+
  # geom_vline(xintercept = 0.1)+
  scale_y_continuous(limits = c(-0.05,1.05),expand=c(0,0), breaks =seq(0,1, by=0.25))+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=30), axis.text = element_text(size=30), legend.title = element_blank(), legend.position = c(0.5,0.9), legend.direction = "horizontal", legend.background = element_rect(fill=NA, color="white"), legend.key = element_rect(fill=NA, color="white"), legend.text = element_text(size = 25))+
  scale_color_manual(breaks=c(0,1,2, 3), values=c("blue3", "purple", "#009E73", "brown2"), labels=c("single-copy pTP  ","single-copy pFP  ", "multicopy pFP  ", "multicopy pTP"))+
  scale_fill_manual(breaks=c(0,1,2, 3), values=c("blue3", "purple", "#009E73", "brown2"), labels=c("single-copy pTP  ","single-copy pFP  ", "multicopy pFP  ", "multicopy pTP"))+
  scale_shape_manual(breaks=c(0,1,2,3), values = c(0,1,2,3), labels=c("single-copy pTP  ","single-copy pFP  ", "multicopy pFP  ", "multicopy pTP"))+
  geom_segment(x=-0.05, xend=-0.05, y=0, yend=1, color="black")+
  geom_segment(x=0, xend=1, y=-0.05, yend=-0.05, color="black")+
  guides(color = guide_legend(override.aes = list(size=10, alpha=1)))
