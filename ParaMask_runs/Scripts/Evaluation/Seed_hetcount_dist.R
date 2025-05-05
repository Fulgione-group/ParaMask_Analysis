
# Comparison of number of heterozygotes at seeds and at nearest seeds of false negatives in Arabis (single-copy SNPs overlapping SVs)

clusters_c1<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr1.clusters.txt", header = T)
clusters_c2<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr2.clusters.txt", header = T)
clusters_c3<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr3.clusters.txt", header = T)
clusters_c4<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr4.clusters.txt", header = T)
clusters_c5<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr5.clusters.txt", header = T)
clusters_c6<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr6.clusters.txt", header = T)
clusters_c7<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr7.clusters.txt", header = T)
clusters_c8<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.chr8.clusters.txt", header = T)

head(clusters_c1)
clusters<- rbind(clusters_c1, clusters_c2, clusters_c3, clusters_c4, clusters_c5, clusters_c6, clusters_c7, clusters_c8)
head(clusters)
het<-read.table("ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES0304_run_finalEMresults.het", header = T)
head(clusters)
head(het)

s<-clusters$hetGenOfLastSeed[clusters$clusterCause=="seed"]
l<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l[6]
l<- as.vector(l)
min(l)
x <- "1:2"
strsplit(x, split = ":")[[1]]
df_as <- data.frame(value = l, type= "All_Seeds")

# Create histogram with 86 bins
ggplot(df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density Plot",
       x = "Value", y = "Density") +
  theme_minimal()

ggplot(df, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Relative Frequencies",
       x = "Value", y = "Relative Frequency") +
  theme_minimal()

### Number of heterozygotes at nearest seeds of false negatives in Arabis (single-copy SNPs overlapping SVs of ES03)

fn_ES03 <- read.table("/ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES03_FN_SNPs_hetGen_of_last_seed_SNPs.txt", header = T)
head(fn_ES03)

s<-fn_ES03$hetGenOfLastSeed[fn_ES03$clusterCause=="seed"]
l_ES03<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l_ES03[6]
l_ES03<- as.vector(l_ES03)
min(l_ES03)
df_ES03 <- data.frame(value = l_ES03, type= "FN_ES03")

# Create histogram with 86 bins
ggplot(df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density Plot",
       x = "Value", y = "Density") +
  theme_minimal()

ggplot(df, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Relative Frequencies",
       x = "Value", y = "Relative Frequency") +
  theme_minimal()

### Number of heterozygotes at nearest seeds of false negatives in Arabis (single-copy SNPs overlapping SVs of ES04)

fn_ES04 <- read.table("/ParaMask_additional_data/ParaMask_results/SpanishArabisAlpina/ES04_FN_SNPs_hetGen_of_last_seed_SNPs.txt", header = T)
head(fn_ES04)

s<-fn_ES04$hetGenOfLastSeed[fn_ES04$clusterCause=="seed"]
l_ES04<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l_ES04[6]
l_ES04<- as.vector(l_ES04)
min(l_ES04)
df_ES04 <- data.frame(value = l_ES04, type= "FN_ES04")

df_all <- rbind(df_as, df_ES03, df_ES04)
# Create histogram with 86 bins
ggplot(df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density Plot",
       x = "Value", y = "Density") +
  theme_minimal()

ggplot(df, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Relative Frequencies",
       x = "Value", y = "Relative Frequency") +
  theme_minimal()

mean(c(df_ES04$value,df_ES03$value))
median(c(df_ES04$value,df_ES03$value))

ggplot(df_all, aes(x = value, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, position = "identity", alpha = 0.5, color = "white") +
  labs(title = "Overlayed Relative Frequency Histograms",
       x = "Value", y = "Relative Frequency") +
  theme_minimal()

max(l)
alpina_plot<- ggplot(df_all, aes(x = value, color =type, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(x = "Number of heterozygotes at seeds", y = "Density") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), limits = c(0,0.3))+
  scale_x_continuous(limits = c(1,85), breaks = c(1,30,60,85))+
  geom_segment(y=-0.015,yend=-0.015, x=1, xend=85, color="black")+
  geom_segment(x=-3.1,xend=-3.1, y=0, yend=0.3, color="black")+
  scale_fill_manual(breaks=c("All_Seeds", "FN_ES03", "FN_ES04"), labels=c("all seeds", "false negative SVs ES03", "false negative SVs ES04"), values =c("brown2", "#1F78B4", "#A6CEE3")
 )+
  scale_color_manual(breaks=c("All_Seeds", "FN_ES03", "FN_ES04"), labels=c("all seeds", "false negative SVs ES03", "false negative SVs ES04"), values =c("brown2", "#1F78B4", "#A6CEE3")
  )+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, color = "white"), axis.title=element_text(size=40), legend.title = element_blank(), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5), legend.position = c(0.6,0.7))


# Species comparison
###load salmon and sinapis data
clusters_ps<-c()
i<-1
for(i in 1:26){
  print(i)
  tmp<-read.table(paste0("ParaMask_additional_data/ParaMask_results/PinkSalmon/PinkSalmon_nosex_EM2.7_EMresults_chr", as.character(i),".clusters.txt"), header = T)
  clusters_ps <- rbind(clusters_ps, tmp)
}


s<-clusters_ps$hetGenOfLastSeed[clusters_ps$clusterCause=="seed"]
l<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l[6]
l<- as.vector(l)
min(l)
df_ps <- data.frame(value = l, type= "O. gorbuscha")


clusters_ls<-c()

for(i in 1:28){
  print(i)
  tmp<-read.table(paste0("ParaMask_additional_data/ParaMask_results/WoodWhiteButterfly/WB_EM_2.7.1_EMresults.chr", as.character(i),".clusters.txt"), header = T)
  clusters_ls <- rbind(clusters_ls, tmp)
}

s<-clusters_ls$hetGenOfLastSeed[clusters_ls$clusterCause=="seed"]
l<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l[6]
l<- as.vector(l)
min(l)
df_ls <- data.frame(value = l, type= "L. sinapis")

mean(df_as$value)
median(df_as$value)

hist(df_ls$value)
df_allspecies <- rbind(df_as, df_ps, df_ls)
max(df_ls$value)
ggplot(df_ls, aes(x = value, color =type, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(x = "Number of heterozygotes at seeds", y = "Density") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), limits = c(0,0.3))+
  scale_x_continuous(limits = c(1,78), breaks = c(1,30,60,85))+
  geom_segment(y=-0.015,yend=-0.015, x=1, xend=85, color="black")+
  geom_segment(x=-3.1,xend=-3.1, y=0, yend=0.3, color="black")+
  scale_fill_manual(breaks=c("All_Seeds", "O. gorbuscha", "L. sinapis"), labels=c(expression(italic("A. alpina")), expression(italic("O. gorbuscha")), expression(italic("L. sinapis"))), values =c("brown2", "#009E73", "#0072B2"))+
  scale_color_manual(breaks=c("All_Seeds", "O. gorbuscha", "L. sinapis"), labels=c(expression(italic("A. alpina")), expression(italic("O. gorbuscha")), expression(italic("L. sinapis"))), values =c("brown2", "#009E73", "#0072B2"))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, color = "white"), axis.title=element_text(size=40), legend.title = element_blank(), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5), legend.position = c(0.6,0.7))

species_plot<-ggplot(df_allspecies, aes(x = value, color =type, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(x = "Number of heterozygotes at seeds", y = "Density") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), limits = c(0,0.3))+
  scale_x_continuous(limits = c(1,85), breaks = c(1,30,60,85))+
  geom_segment(y=-0.015,yend=-0.015, x=1, xend=85, color="black")+
  geom_segment(x=-3.1,xend=-3.1, y=0, yend=0.3, color="black")+
  scale_fill_manual(breaks=c("All_Seeds", "O. gorbuscha", "L. sinapis"), labels=c(expression(italic("A. alpina")), expression(italic("O. gorbuscha")), expression(italic("L. sinapis"))), values =c("brown2", "#009E73", "#0072B2"))+
  scale_color_manual(breaks=c("All_Seeds", "O. gorbuscha", "L. sinapis"), labels=c(expression(italic("A. alpina")), expression(italic("O. gorbuscha")), expression(italic("L. sinapis"))), values =c("brown2", "#009E73", "#0072B2"))+
  theme(axis.text.x= element_text(size=40),axis.text.y = element_text(size=40), panel.background = element_rect(fill=NA, color = "white"), axis.title=element_text(size=40), legend.title = element_blank(), legend.key.size = unit(2, "cm"),legend.key = element_rect(fill="white"), legend.text = element_text(size = 40), axis.ticks.length.y  =unit(0.5, "cm"), axis.ticks.length.x  =unit(0.5, "cm"), plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size=40, hjust = 0.5, vjust = 1.5), legend.position = c(0.6,0.7))



library(cowplot)
plot_grid(plotlist = list(alpina_plot, species_plot), ncol = 2, rel_widths = c(0.5,0.5))


sum(l==1 | l==2 | l==3)/length(l)



sum(l==1 | l==2 | l==3)/length(l)

clusters[clusters$clusterCause=="seed",][l==1,]
pos <- clusters$position[clusters$clusterCause=="seed"][l==1]
pos <- clusters$position[clusters$clusterCause=="seed"][l==2]
pos <- clusters$position[clusters$clusterCause=="seed"][l==3]
pos <- clusters$position[clusters$clusterCause=="seed"][l==4]

het[het$Chromosome=="chr1" & het$Position %in%pos,]
test <-het[het$Chromosome=="chr1" & het$Position %in%pos,]


s<-clusters$hetGenOfLastSeed[clusters$clusterCause=="bridge"]
l<- sapply(s, FUN = function(x){length(strsplit(x, split = ":")[[1]])})
l[6]
l<- as.vector(l)
min(l)
x <- "1:2"
strsplit(x, split = ":")[[1]]
df <- data.frame(value = l)

# Create histogram with 86 bins
ggplot(df, aes(x = value)) +
  geom_histogram(binwidth  = 1, fill = "steelblue", color = "white") +
  labs(title = "Histogram with 86 Bins",
       x = "Value", y = "Count") +
  theme_minimal()

