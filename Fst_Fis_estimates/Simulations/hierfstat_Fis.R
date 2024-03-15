

# install.packages("hierfstat")
library(hierfstat)
library(vcfR)
library(adegenet)
library(ggplot2)



vcf<-read.vcfR("~/Data/Paramask/Sim_final/HW_5PercentParalogs/VCFs/Simulations_SV_SC_rep1_mat_cpos.pruned.plink.vcf")
geneind <- vcfR2genind(x = vcf, return.alleles = T)
pop(geneind) <- as.factor(rep("tmp",100))
b<-basic.stats(geneind)
head(b$perloc)
Ho(geneind)
Hs(geneind)

vcf<-read.vcfR("~/Data/Paramask/Sim_final/HW_0PercentParalogs/VCFs/Simulations_SV_SC_rep1_mat_cpos.pruned.plink.p1.vcf")
geneind <- vcfR2genind(x = vcf, return.alleles = T)
pop(geneind) <- as.factor(rep("tmp",100))
Ho(geneind)
Hs(geneind)

vcf<-read.vcfR("~/Data/Paramask/Sim_final/HW_0PercentParalogs/VCFs/Simulations_SV_SC_rep1_mat_cpos.pruned.plink.p2.vcf")
geneind <- vcfR2genind(x = vcf, return.alleles = T)
pop(geneind) <- as.factor(rep("tmp",100))
Ho(geneind)
Hs(geneind)

##

Fis_tab<-c()
i<-1
for(i in 1:3){
  vcf<-read.vcfR(paste("~/Data/Paramask/Sim_final/HW_0PercentParalogs/VCFs/Simulations_SV_SC_rep", as.character(i),"_mat_cpos.pruned.plink.vcf", sep=""))
  geneind <- vcfR2genind(x = vcf, return.alleles = T)
  pop(geneind) <- as.factor(rep("tmp",100))
  bs<-basic.stats(geneind)
  f <- bs$overall[9]
  Fis_tab<-rbind(Fis_tab, c(0,"SC", i, f))
}


for(p in seq(5,50, by=5)){
  for(t in c("","_sc", "_ParaOut")){
    for(i in 1:3){
      vcf<-read.vcfR(paste("~/Data/Paramask/Sim_final/HW_",as.character(p),"PercentParalogs/VCFs/Simulations_SV_SC_rep", as.character(i),"_mat_cpos", t,".pruned.plink.vcf", sep=""))
      geneind <- vcfR2genind(x = vcf, return.alleles = T)
      pop(geneind) <- as.factor(rep("tmp",100))
      bs<-basic.stats(geneind)
      f <- bs$overall[9]
      Fis_tab<-rbind(Fis_tab, c(p,t, i, f))
    }
  }
}
Fis_tab <- as.data.frame(Fis_tab)
colnames(Fis_tab)<- c("proportion", "type", "replicate", "Fis")
Fis_tab$proportion[-c(1:3)]<- sort(rep(seq(5,50, by=5),9))
Fis_tab$type[-c(1:3)]<- rep(c(rep("collapsed", 3), rep("sc",3),rep("ParaOut",3)),10)
# Fis_tab$type[c(1:3)]<- rep("sc",3)
# write.table(Fis_tab, file = "~/Data/Paramask/Sim_final/Fis_tab.txt", col.names = T, sep = "\t")
Fis_tab<- read.table(file = "~/Data/Paramask/Sim_final/Fis_tab.txt", header = T, sep = "\t")
Fis_tab$Fis<- as.numeric(Fis_tab$Fis)
Fis_tab$proportion<- as.numeric(Fis_tab$proportion)
min(Fis_tab$Fis)
Fis_tab$type
Fis_tab$replicate<- as.factor(Fis_tab$replicate)

ggplot(data = Fis_tab, aes(x = proportion, y=Fis, color=type, shape=replicate))+
  geom_point(size=4,stroke=2)+
  scale_x_continuous(breaks = seq(0,50, by=5), labels = paste(seq(0,50, by=5), "%", sep = ""), limits = c(-2.5,52.5), expand=c(0,0))+
  labs(y= expression(F[IS]), x="Proportion of duplications")+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=16), axis.text = element_text(size=14), legend.title = element_blank(), legend.text = element_text(size=14))+
  scale_y_continuous(breaks = seq( 0.1, -0.6, by=-0.1), limits = c(-0.65,0.15), expand=c(0,0))+
  geom_segment(x=-2.5, xend=-2.5, y=-0.6, yend=0.1, color="black")+
  geom_segment(x=0, xend=50, y=-0.65, yend=-0.65, color="black")+
  scale_color_brewer(palette = "Dark2",breaks=c("collapsed", "ParaOut", "sc"), labels=c("collapsed", "masked", "simulated single copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))

seq(-0.6, 0.1, by=0.1)
seq( 0.1, -0.6, by=-0.1)


  
plot_F<-ggplot(data = Fis_tab, aes(x = proportion, y=Fis, color=type, shape=replicate))+
  geom_point(size=6,stroke=2)+
  scale_x_continuous(breaks = seq(0,50, by=5), limits = c(-4,54), expand=c(0,0))+
  labs(y= expression(F[IS]), x="Proportion of duplications (%)")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x = element_blank(), axis.title.y = element_text(size=40, vjust = 1), axis.text = element_text(size=40), legend.title = element_blank(), legend.text = element_text(size=14), legend.position = "", plot.margin = margin(t=90, r=1, l=1, b=1))+
  scale_y_continuous(breaks = seq( 0, -0.5, by=-0.25), limits = c(-0.65,0.15), expand=c(0,0))+
  geom_segment(x=-4, xend=-4, y=-0.5, yend=0, color="black")+
  geom_segment(x=0, xend=50, y=-0.65, yend=-0.65, color="black")+
  scale_color_manual(values=c("#920294", "#7eb302", "blue3"),breaks=c("collapsed", "ParaOut", "sc"), labels=c("collapsed", "masked", "simulated single copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))

#4by4
plot_F<-ggplot(data = Fis_tab, aes(x = proportion, y=Fis, color=type, shape=replicate))+
  geom_point(size=6,stroke=2)+
  scale_x_continuous(breaks = seq(0,50, by=5), limits = c(-4,54), expand=c(0,0))+
  labs(y= expression(F[IS]), x="Proportion of duplications (%)")+
  theme(axis.ticks.length=unit(.25, "cm"),panel.background = element_rect(fill=NA, color="white"), axis.title.x = element_text(size=40, vjust = -0.75), axis.title.y = element_text(size=40, vjust = 1), axis.text = element_text(size=40), legend.title = element_blank(), legend.text = element_text(size=14), legend.position = "", plot.margin = margin(t=90, r=1, l=60, b=60))+
  scale_y_continuous(breaks = seq( 0, -0.5, by=-0.25), limits = c(-0.65,0.15), expand=c(0,0))+
  geom_segment(x=-4, xend=-4, y=-0.5, yend=0, color="black")+
  geom_segment(x=0, xend=50, y=-0.65, yend=-0.65, color="black")+
  scale_color_manual(values=c("#920294", "#7eb302", "blue3"),breaks=c("collapsed", "ParaOut", "sc"), labels=c("collapsed", "masked", "simulated single copy"))+
  scale_shape_manual(values = c(0,1,2),labels=c("replicate 1", "replicate 2", "replicate 3"))


Fis_bias <- c()
Fis_corrected_bias <- c()
i <-5
for(i in seq(5,50, by=5)){
  SC<- Fis_tab[Fis_tab$proportion==i & Fis_tab$type=="sc",]
  WG <- Fis_tab[Fis_tab$proportion==i & Fis_tab$type=="collapsed",]
  PO <- Fis_tab[Fis_tab$proportion==i & Fis_tab$type=="ParaOut",]
  Fis_bias <- rbind(Fis_bias,c(i, mean(WG$Fis-SC$Fis)))
  Fis_corrected_bias <- rbind(Fis_corrected_bias,c(i, mean(PO$Fis-SC$Fis)))
}
Fis_bias
Fis_corrected_bias


