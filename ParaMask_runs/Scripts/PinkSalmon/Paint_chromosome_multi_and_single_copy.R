


bed<-read.table(file = "/PinkSalmon_nosex_EM2.7.1_EMresults_chrall.finalClass.bed", sep = "\t",  header=T)

head(bed)


bed$type.0.single.copy.1.multi.copy<- factor(bed$type.0.single.copy.1.multi.copy, levels = c(0,1))
bed$Chromosome <- factor(bed$Chromosome, levels = unique(bed$Chromosome))
max(bed$End)
head(bed$Chromosome)

ggplot(data = bed, aes(group=Chromosome))+
  geom_segment( mapping = aes(x = Start, xend=End, y=Chromosome, yend=Chromosome,color=type.0.single.copy.1.multi.copy), size=4)+
  scale_x_continuous(limits = c(0, 122000000), expand = c(0,2000000), breaks=seq(0, 120000000, by=10000000), labels=seq(0,120, by=10))+
  scale_y_discrete(expand = c(0,1))+
  labs(x="Position [Mbp]", y="Chromosome")+
  scale_color_manual(breaks=c(0,1), values = c("blue3","brown2"), labels=c("single-copy    ","multicopy"))+
  geom_segment(x=0, xend=120000000, y= 0, yend=0, position = position_nudge(x=0,y = 0))+
  theme(axis.title.y=element_text(size=25, vjust = 3),axis.title.x=element_text(size=25, vjust = -1) ,panel.background = element_rect(fill=NA, color="white"), legend.position = "top",aspect.ratio = 1, axis.ticks.y = element_blank(), axis.text = element_text(size=20),legend.key.width = unit(2, 'cm'), legend.key.height  = unit(1, 'cm'), legend.key = element_rect(fill = "white"), legend.text = element_text(size=20), plot.margin = margin(t = 10, r = 5, b = 15, l = 15))+
  labs(x="Position [Mbp]", y="Chromosome", color="")
  
