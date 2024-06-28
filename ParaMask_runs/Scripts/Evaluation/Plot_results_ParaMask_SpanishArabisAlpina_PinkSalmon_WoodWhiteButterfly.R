library(ggplot2)


bed<-read.table(file = "ES0304_run_finalEMresults.allChr.finalClass.multicopy.bed", sep = "\t", header = F)
bed<-read.table(file = "PinkSalmon_nosex_EM2.7_EMresults_chrall.finalClass.multicopy.bed", sep = "\t", header = F)
bed<-read.table(file = "WB_EM_2.7.1_EMresults.chrall.finalClass.multicopy.bed", sep = "\t", header = F)



colnames(bed)<- c("chr", "start", "end", "type", "nSNPs", "cluster")
head(bed)

bed$length <- bed$end- bed$start + 1

sum(bed$length)
max(bed$length)
mean(bed$length)
sum(bed$length)/311583516

mean(bed$nSNPs)

bed$length2 <- bed$length
bed$length2[bed$length2>=3000]<-3000
bed$length2[bed$length2>=20000]<-20000


max(bed$length2)
min(bed$length2)
sum(is.na(bed$length2))
sum(bed$length2==2)

bed$length2 <- as.numeric(bed$length2)

bed$length3 <- bed$length2-10  
bed$length3 <- bed$length2-20

hist(bed$length)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
tabulate(bed$length)
Mode(bed$length)
Mode(bed$length[bed$length<550])
Mode(bed$length[bed$length<301])
Mode(bed$length[bed$length<159])

Mode(bed$length[bed$length>360])

hist(bed$length3)

ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 20, color="black", fill="brown2", linewidth=0.2)+
  scale_x_continuous(limits = c(-10,3000), expand=c(0,90), breaks=c(-10,980,1980, 2980), labels=c(0,1,2,expression("">=3)))+
  scale_y_continuous(expand=c(0,0), limits = c(-0.005,0.08), breaks=c(0,0.04,0.08), labels=c(0,4,8))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  geom_segment(x=-10, xend=2980, y=-0.005, yend=-0.005)+
  geom_segment(x=-100, xend=-100, y=0, yend=0.08)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_text(size=40, vjust = -0.5, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 4) , axis.text.y  = element_text(size=40), axis.text.x  = element_text(size=40, vjust=-0.25), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))


##inlet for Pink Salmon

ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 20, color="black", fill="brown2", linewidth=0.2)+
  scale_x_continuous(limits = c(-10,1000), expand=c(0,0), breaks=c(-10,980,1980, 2980), labels=c(0,1,2,expression("">=3)))+
  scale_y_continuous(expand=c(0,0), limits = c(0,0.10), breaks=c(0,0.04,0.08), labels=c(0,4,8))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  # geom_segment(x=-10, xend=2980, y=-0.005, yend=-0.005)+
  # geom_segment(x=-100, xend=-100, y=0, yend=0.08)+
  theme(axis.ticks=element_blank(),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_blank(),  axis.title.y = element_blank() , axis.text.y  = element_blank(), axis.text.x  = element_blank(), legend.title = element_blank(), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))
##

ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 20, color="black", fill="brown2", linewidth=0.2)+
  scale_x_continuous(limits = c(-10,3000), expand=c(0,90), breaks=c(-10,980,1980, 2980), labels=c(0,1,2,expression("">=3)))+
  scale_y_continuous(expand=c(0,0), limits = c(-0.005,0.5), breaks=seq(0,0.4, by=0.2), labels=seq(0,40, by=20))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  geom_segment(x=-10, xend=2980, y=-0.005, yend=-0.005)+
  geom_segment(x=-100, xend=-100, y=0, yend=0.4)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_text(size=40, vjust = -0.5, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 4) , axis.text.y  = element_text(size=40), axis.text.x  = element_text(size=40, vjust=-0.25), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))

#wb
ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 20, color="black", fill="brown2", linewidth=0.2)+
  scale_x_continuous(limits = c(-10,3000), expand=c(0,90), breaks=c(-10,980,1980, 2980), labels=c(0,1,2,expression("">=3)))+
  scale_y_continuous(expand=c(0,0), limits = c(-0.005,0.1), breaks=seq(0,0.1, by=0.05), labels=seq(0,10, by=5))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  geom_segment(x=-10, xend=2980, y=-0.005, yend=-0.005)+
  geom_segment(x=-100, xend=-100, y=0, yend=0.4)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_text(size=40, vjust = -0.5, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 4) , axis.text.y  = element_text(size=40), axis.text.x  = element_text(size=40, vjust=-0.25), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))

##20kbp cutoff Salmon
min(bed$length3)
ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth = 20, color="black", fill="brown2", linewidth=0.05)+
  scale_x_continuous(limits = c(-1000,21000), expand=c(0,90), breaks=c(-10,10980,19980), labels=c(0,10,expression("">=20)))+
  scale_y_continuous(expand=c(0,0), limits = c(-0.005,0.1), breaks=c(0,0.05,0.1), labels=c(0,5,10))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  geom_segment(x=-10, xend=19980, y=-0.005, yend=-0.005)+
  geom_segment(x=-1000, xend=-1000, y=0, yend=0.4)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_text(size=40, vjust = -0.5, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 4) , axis.text.y  = element_text(size=40), axis.text.x  = element_text(size=40, vjust=-0.25), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))

ggplot(data=bed, aes(x=length3))+
  geom_histogram(aes(y = after_stat(count / sum(count))),binwidth =40, color="black", fill="brown2", linewidth=0.05)+
  scale_x_continuous(limits = c(-20,20000), expand=c(0,990), breaks=c(-10,10980,19980), labels=c(0,10,expression("">=20)))+
  scale_y_continuous(expand=c(0,0), limits = c(-0.005,0.1), breaks=c(0,0.05,0.1), labels=c(0,5,10))+
  labs(y=expression(paste("Frequency (" ,10^{-2},")")), x= "Length (kbp)")+
  geom_segment(x=-10, xend=19980, y=-0.005, yend=-0.005)+
  geom_segment(x=-1000, xend=-1000, y=0, yend=0.4)+
  geom_segment(x=0, xend=0, y=0, yend=-0.0005)+
  geom_segment(x=1000, xend=1000, y=0, yend=-0.0005)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x  =  element_text(size=40, vjust = -0.5, hjust=0.5),  axis.title.y = element_text(size=40, vjust = 4) , axis.text.y  = element_text(size=40), axis.text.x  = element_text(size=40, vjust=-0.25), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))


het <- read.table(file = "ES0304_run_finalEMresults.chrall.finalClass.het", sep = "\t", header = T)
het <- read.table(file = "~/Data/Paramask/Validation_ParaMask/PinkSalmon/PinkSalmon_nosex_EM2.7_EMresults_chrall.finalClass.het", sep = "\t", header = T)
het <- read.table(file = "~/Data/Paramask/Validation_ParaMask/WB/WB_EM_2.7.1_EMresults.chrall.finalClass.het", sep = "\t", header = T)

max(het$Non.missing)
head(het)
sum(het$finalClass==1)/nrow(het)

EM_par_seeds<-sum(het$EM_class==2 & het$allele.deviation.seed==0)
EM_sc_seeds<-sum(het$EM_class==0)
AR_par_seeds <- sum(het$allele.deviation.seed==1)

par_SNPs<-sum(het$finalClass==1)
sc_SNPs<-sum(het$finalClass==0)
total_snps <- par_SNPs + sc_SNPs

EM_par_seeds/par_SNPs
AR_par_seeds/par_SNPs



rd_results <- data.frame(N=c(par_SNPs, EM_par_seeds, AR_par_seeds, sc_SNPs, EM_sc_seeds), type=c("par", "par", "par", "sc", "sc"), classifier=c("overall", "EM", "AR", "overall", "EM"))

rd_results$Freq <- rd_results$N/total_snps

rd_results$plotfac <- paste(rd_results$type,rd_results$classifier, sep="_")
rd_results$plotfac <- factor(rd_results$plotfac, levels=rd_results$plotfac)

ggplot(data = rd_results, aes(x = plotfac, y=Freq))+
  labs(x="Proportion of duplications (%)", y="Proportion of SNPs")+
  geom_bar(aes(fill=plotfac , alpha=plotfac),stat="identity", position = "dodge", width = 0.5,color="black", size=1.1)+
  # geom_text(data = text_labels, aes(x = x, y = Correct, label =classifier),size=5, angle=45)+
  # geom_point(data = bar_labels, aes(y=Correct, x=x, color=type, shape=Replicate), size=6, stroke=2)+
  scale_x_discrete(limits=c(levels(rd_results$plotfac))[5:1], breaks=levels(rd_results$plotfac), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_y_continuous(position = "right",breaks=seq(0,0.8, by=0.4), limits=c(-0.025,0.825), expand=c(0,0))+
  scale_fill_manual(values=c("brown2", "brown2", "brown2", "blue3","blue3"),breaks=c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_alpha_manual(breaks = c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), values = c(1, 0.6, 0.2, 1, 0.6),labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  geom_segment(data = rd_results,aes(x="par_overall", xend="par_overall", y=0, yend=0.8), color="black", position = position_nudge(x = 0.6))+
  geom_segment(data = rd_results, aes(x="sc_EM", xend="par_overall", y=-0.025, yend=-0.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.2, yend=1.2), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.2, yend=1.2), color="black")+
  # annotate(geom = "text",x = "0.1_0", y = 0,label="F[IS] == 0", parse=T,vjust =3.5, hjust=0, size=12, angle =90)+
  # annotate(geom = "text",x = "0.1_0.9", y = 0,label="F[IS] == 0.9", parse=T,vjust =3.5, size=12, angle=90)+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x.top  =  element_text(size=40, vjust = 4, hjust=0.5),  axis.title.y = element_blank(), axis.text.y = element_text(size=40), axis.text.x = element_text(size=40, vjust = 0.5), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))+
  coord_flip()+
  guides(
    fill = guide_legend(override.aes = list( linetype = 0), nrow=3, byrow = T)
  )
  

ggplot(data = rd_results, aes(x = plotfac, y=Freq))+
  labs(x="Proportion of duplications (%)", y="Proportion of SNPs")+
  geom_bar(aes(fill=plotfac , alpha=plotfac),stat="identity", position = "dodge", width = 0.5,color="black", size=1.1)+
  # geom_text(data = text_labels, aes(x = x, y = Correct, label =classifier),size=5, angle=45)+
  # geom_point(data = bar_labels, aes(y=Correct, x=x, color=type, shape=Replicate), size=6, stroke=2)+
  scale_x_discrete(limits=c(levels(rd_results$plotfac))[5:1], breaks=levels(rd_results$plotfac), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"), expand=c(0,0))+
  scale_y_continuous(position = "right",breaks=seq(0,1, by=0.5), limits=c(0,1.03), expand=c(0,0))+
  scale_fill_manual(values=c("brown2", "brown2", "brown2", "blue3","blue3"),breaks=c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_alpha_manual(breaks = c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), values = c(1, 0.6, 0.2, 1, 0.6),labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  geom_segment(data = rd_results,aes(x="par_overall", xend="par_overall", y=0, yend=1), color="black", position = position_nudge(x = 0.6))+
  # geom_segment(data = rd_results, aes(x="sc_EM", xend="par_overall", y=-0.025, yend=-0.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.2, yend=1.2), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.2, yend=1.2), color="black")+
  annotate(geom = "text",x = "par_overall", y = rd_results$Freq[1]+0.02,label="mc: overall", parse=F, hjust=0, size=12, angle =0)+
  annotate(geom = "text",x = "par_EM", y = rd_results$Freq[2]+0.02,label="mc: EM", parse=F, hjust=0, size=12, angle =0)+
  annotate(geom = "text",x = "par_AR", y = rd_results$Freq[3]+0.02,label="mc: RRD", parse=F, hjust=0, size=12, angle =0)+
  # annotate(geom = "text",x = "sc_overall", y = rd_results$Freq[4]+0.02,label="sc: overall", parse=F, hjust=0, size=12, angle =0)+
  annotate(geom = "text",x = "sc_EM", y = rd_results$Freq[5]+0.02,label="sc: EM", parse=F, hjust=0, size=12, angle =0)+
  
  # annotate(geom = "text",x = "0.1_0.9", y = 0,label="F[IS] == 0.9", parse=T,vjust =3.5, size=12, angle=90)+
  theme(axis.ticks.length=unit(.25, "cm"), axis.ticks.y = element_blank(),panel.background = element_rect(fill=NA, color="white"), axis.title.x.top  =  element_text(size=36, vjust = 2, hjust=0.5),  axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=40, vjust = 0.5), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))+
  coord_flip()+
  guides(
    fill = guide_legend(override.aes = list( linetype = 0), nrow=3, byrow = T)
  )

ggplot(data = rd_results, aes(x = plotfac, y=Freq))+
  labs(x="Proportion of duplications (%)", y="Proportion of SNPs")+
  geom_bar(aes(fill=plotfac , alpha=plotfac),stat="identity", position = "dodge", width = 0.5,color="black", size=1.1)+
  # geom_text(data = text_labels, aes(x = x, y = Correct, label =classifier),size=5, angle=45)+
  # geom_point(data = bar_labels, aes(y=Correct, x=x, color=type, shape=Replicate), size=6, stroke=2)+
  scale_x_discrete(limits=c(levels(rd_results$plotfac))[5:1], breaks=levels(rd_results$plotfac), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"), expand=c(0,0))+
  scale_y_continuous(position = "right",breaks=seq(0,1, by=0.5), limits=c(0,1.03), expand=c(0,0))+
  scale_fill_manual(values=c("brown2", "brown2", "brown2", "blue3","blue3"),breaks=c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  scale_alpha_manual(breaks = c("par_overall","par_EM", "par_AR","sc_overall","sc_EM" ), values = c(1, 0.6, 0.2, 1, 0.6),labels=c("multicopy: overall","multicopy: EM step", "multicopy: RRD step", "single-copy: overall","single-copy: EM step"))+
  geom_segment(data = rd_results,aes(x="par_overall", xend="par_overall", y=0, yend=1), color="black", position = position_nudge(x = 0.6))+
  # geom_segment(data = rd_results, aes(x="sc_EM", xend="par_overall", y=-0.025, yend=-0.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.025, yend=1.025), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0.9", xend="0.1_0.9", y=1.2, yend=1.2), color="black")+
  # geom_segment(data = stab, aes(x="0.5_0", xend="0.1_0", y=1.2, yend=1.2), color="black")+
  annotate(geom = "text",x = "par_overall", y = rd_results$Freq[1]+0.02,label="mc: overall", parse=F, hjust=0, size=10, angle =0)+
  annotate(geom = "text",x = "par_EM", y = rd_results$Freq[2]+0.02,label="mc: EM", parse=F, hjust=0, size=10, angle =0)+
  annotate(geom = "text",x = "par_AR", y = rd_results$Freq[3]+0.02,label="mc: RRD", parse=F, hjust=0, size=10, angle =0)+
  # annotate(geom = "text",x = "sc_overall", y = rd_results$Freq[4]+0.02,label="sc: overall", parse=F, hjust=0.3, size=10, angle =0)+
  annotate(geom = "text",x = "sc_EM", y = rd_results$Freq[5]+0.02,label="sc: EM", parse=F, hjust=0, size=10, angle =0)+
  
  # annotate(geom = "text",x = "0.1_0.9", y = 0,label="F[IS] == 0.9", parse=T,vjust =3.5, size=12, angle=90)+
  theme(axis.ticks.length=unit(.25, "cm"), axis.ticks.y = element_blank(),panel.background = element_rect(fill=NA, color="white"), axis.title.x.top  =  element_text(size=36, vjust = 2, hjust=0.5),  axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=40, vjust = 0.5), legend.title = element_blank(), legend.text = element_text(size=30), legend.position = "", plot.margin = margin(t=30, r=30, l=60, b=30), legend.key.size = unit(1.7, "lines"),legend.spacing.x = unit(0.5, "cm"), legend.key.width = unit(2, "cm"), legend.key = element_rect(color="black", size=1.1), legend.spacing.y = unit(0.3, "cm"))+
  coord_flip()+
  guides(
    fill = guide_legend(override.aes = list( linetype = 0), nrow=3, byrow = T)
  )

#


