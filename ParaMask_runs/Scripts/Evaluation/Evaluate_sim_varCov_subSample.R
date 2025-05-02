colSdApply <- function(x, ...)apply(X=x, MARGIN=2, FUN=sd, ...)

sim_stats <- c()
sim_stats_av <- c()
i<- 1
j<-3
k<-10
for(k in c(15,20,50,100)){
  for(j in c(3,5,8,12,20)){
    sim_stats_tmp <- c()
    for(i in 1:3){
      paraOut<-read.table(file = paste0("~/Data/Paramask/Sim_final/Sim_varCov_SubSample/HW_10PercentParalogs/Sim0.1HW_cov",as.character(j),"_sub",as.character(k),"_rep",as.character(i),"_EMresults.finalClass.het"), header=T)
      SVtab <- read.table(file = paste0("~/Data/Paramask/Test_sim_rCNV/HW_prop/rep",as.character(i),"/Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
      colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")
      paraOut$Position<- as.numeric(paraOut$Position)
      SVtab$POS<- as.numeric(SVtab$POS)
      SVtab$para_call <- "U"
      SVtab$para_call[SVtab$POS %in% paraOut$Position] <-ifelse(test = paraOut$finalClass==0, yes = "SC", no = "SV")
      SVtab<- SVtab[SVtab$para_call!="U",]
      a<- sum(SVtab$true_cn=="SC" &SVtab$para_call=="SC")/sum(SVtab$true_cn=="SC")
      b<- sum(SVtab$true_cn=="SC" &SVtab$para_call=="SV")/sum(SVtab$true_cn=="SC")
      c<- sum(SVtab$true_cn=="SV" &SVtab$para_call=="SC")/sum(SVtab$true_cn=="SV")
      d<- sum(SVtab$true_cn=="SV" &SVtab$para_call=="SV")/sum(SVtab$true_cn=="SV")
      e<- sum(SVtab$true_cn==SVtab$para_call)/nrow(SVtab)
      sim_stats_tmp <- rbind(sim_stats_tmp, c(a,b,c,d,e,j,k))
    }
    sim_stats <- rbind(sim_stats,sim_stats_tmp)
    sim_stats_av <- rbind(sim_stats_av, c(colMeans(sim_stats_tmp[,c(1:5)]),colSdApply(sim_stats_tmp[,c(1:5)]),j,k))
  }
}
colnames(sim_stats) <- c("Correct_SC", "Incorrect_SC","Incorrect_SV", "Correct_SV", "Total_Recall", "Coverage", "Samplesize")
colnames(sim_stats_av) <- c("M_Correct_SC", "M_Incorrect_SC","M_Incorrect_SV", "M_Correct_SV", "M_Total_Recall","SD_Correct_SC", "SD_Incorrect_SC","SD_Incorrect_SV", "SD_Correct_SV", "SD_Total_Recall", "Coverage", "Samplesize")


# plot_tab<- rbind(cbind(sim_stats_av[,1],sim_stats_av[,6],sim_stats_av[,11], "SC") ,cbind(sim_stats_av[,4],sim_stats_av[,9],sim_stats_av[,11], "SV") , cbind(sim_stats_av[,5],sim_stats_av[,10],sim_stats_av[,11], "Total") )
#
# colnames(plot_tab)<- c("Mean_recall", "SD_recall", "sample_size", "Type")
#
# plot_tab<- as.table(plot_tab)
sim_stats <- as.data.frame(sim_stats)
sim_stats_av <- as.data.frame(sim_stats_av)
# plot_tab<- rbind(cbind(sim_stats[,1],sim_stats[,6], "SC") ,cbind(sim_stats[,4],sim_stats[,6], "SV") , cbind(sim_stats[,5],sim_stats[,6], "Total") )
#
# #
# colnames(plot_tab)<- c("Recall", "Coverage", "Type")
# plot_tab<- as.data.frame(plot_tab)
# plot_tab$Recall <- as.numeric(plot_tab$Recall)
# plot_tab$rep <- as.factor(rep(c(1:3),15))
# plot_tab_HW<- plot_tab

head(sim_stats)
min(sim_stats$Correct_SC)
# plot_tab_HW$Coverage<- factor(plot_tab_HW$Coverage, levels = c("3x","5x","8x","12x", "20x"))

sim_stats_av$Coverage<- factor(sim_stats_av$Coverage, levels= c(3,5,8,12,20))
library(ggplot2)
sim_stats_av_HW <- sim_stats_av
min(sim_stats_av_HW$M_Correct_SC)



HW_sv<-ggplot(sim_stats_av_HW, aes(x = Samplesize, y = M_Correct_SV, color = Coverage, group = Coverage)) +
  geom_point(size = 5, alpha=0.3) +
  geom_line(size=1.1) +
  geom_errorbar(aes(ymin = M_Correct_SV - SD_Correct_SV, ymax = M_Correct_SV +SD_Correct_SV), width = 2, size=1.1) +  # Adjust 'SD' to your column name
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(12.5, 102.5), breaks = c(15, 20, 50, 100)) +
  geom_segment(y = 0, yend = 1, x = 8, xend = 8, color = "black") +
  geom_segment(y = -0.05, yend = -0.05, x = 15, xend = 100, color = "black") +
  labs(x = "", y = "Recall", color = "") +
  scale_color_manual(
    breaks = c(3, 5, 8, 12, 20),
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  theme(
    axis.text.x = element_text(size = 36, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5)
  ) +
  guides(
    shape = "none",
    color = "none"
  )
head(sim_stats_av)
HW_sc<-ggplot(sim_stats_av_HW, aes(x = Samplesize, y = M_Correct_SC, color = Coverage, group = Coverage)) +
  geom_point(size = 5, alpha=0.3) +
  geom_line(size=1.1) +
  geom_errorbar(aes(ymin = M_Correct_SC - SD_Correct_SC, ymax = M_Correct_SC +SD_Correct_SC), width = 2, size=1.1) +  # Adjust 'SD' to your column name
  scale_y_continuous(limits = c(0.985, 1)) +
  scale_x_continuous(limits = c(12.5, 102.5), breaks = c(15, 20, 50, 100)) +
  geom_segment(y = 0.985, yend = 1, x = 8, xend = 8, color = "black") +
  geom_segment(y = 0.98425, yend = 0.98425, x = 15, xend = 100, color = "black") +
  labs(x = "Sample size", y = "Recall", color = "") +
  scale_color_manual(
    breaks = c(3, 5, 8, 12, 20),
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  theme(
    axis.text.x = element_text(size = 36, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5)
  ) +
  guides(
    shape = "none",
    color = "none"
  )

  # Increase legend point size














##############################
#Fis0.9

sim_stats <- c()
sim_stats_av <- c()
i<- 1
j<-3
k<-10
for(k in c(10,15,20,50,100)){
  for(j in c(3,5,8,12,20)){
    sim_stats_tmp <- c()
    for(i in 1:3){
      paraOut<-read.table(file = paste0("~/Data/Paramask/Sim_final/Sim_varCov_SubSample/Fis0.9_10PercentParalogs/Sim0.1Fis0.9_cov",as.character(j),"_sub",as.character(k),"_rep",as.character(i),"_EMresults.finalClass.het"), header=T)
      SVtab <- read.table(file = paste0("~/Data/Paramask/Test_sim_rCNV/selfing/rep",as.character(i),"/Selfing_Simulations_SV_SC_rep",as.character(i),"_mat_cpos.biallelic_true.txt"))
      colnames(SVtab) <- c("Chr", "POS", "ID", "true_cn")
      paraOut$Position<- as.numeric(paraOut$Position)
      SVtab$POS<- as.numeric(SVtab$POS)
      SVtab$para_call <- "U"
      SVtab$para_call[SVtab$POS %in% paraOut$Position] <-ifelse(test = paraOut$finalClass==0, yes = "SC", no = "SV")
      SVtab<- SVtab[SVtab$para_call!="U",]
      a<- sum(SVtab$true_cn=="SC" &SVtab$para_call=="SC")/sum(SVtab$true_cn=="SC")
      b<- sum(SVtab$true_cn=="SC" &SVtab$para_call=="SV")/sum(SVtab$true_cn=="SC")
      c<- sum(SVtab$true_cn=="SV" &SVtab$para_call=="SC")/sum(SVtab$true_cn=="SV")
      d<- sum(SVtab$true_cn=="SV" &SVtab$para_call=="SV")/sum(SVtab$true_cn=="SV")
      e<- sum(SVtab$true_cn==SVtab$para_call)/nrow(SVtab)
      sim_stats_tmp <- rbind(sim_stats_tmp, c(a,b,c,d,e,j,k))
    }
    sim_stats <- rbind(sim_stats,sim_stats_tmp)
    sim_stats_av <- rbind(sim_stats_av, c(colMeans(sim_stats_tmp[,c(1:5)]),colSdApply(sim_stats_tmp[,c(1:5)]),j,k))
  }
}
colnames(sim_stats) <- c("Correct_SC", "Incorrect_SC","Incorrect_SV", "Correct_SV", "Total_Recall", "Coverage", "Samplesize")
colnames(sim_stats_av) <- c("M_Correct_SC", "M_Incorrect_SC","M_Incorrect_SV", "M_Correct_SV", "M_Total_Recall","SD_Correct_SC", "SD_Incorrect_SC","SD_Incorrect_SV", "SD_Correct_SV", "SD_Total_Recall", "Coverage", "Samplesize")


# plot_tab<- rbind(cbind(sim_stats_av[,1],sim_stats_av[,6],sim_stats_av[,11], "SC") ,cbind(sim_stats_av[,4],sim_stats_av[,9],sim_stats_av[,11], "SV") , cbind(sim_stats_av[,5],sim_stats_av[,10],sim_stats_av[,11], "Total") )
#
# colnames(plot_tab)<- c("Mean_recall", "SD_recall", "sample_size", "Type")
#
# plot_tab<- as.table(plot_tab)
sim_stats <- as.table(sim_stats)
sim_stats_av <- as.data.frame(sim_stats_av)
# plot_tab<- rbind(cbind(sim_stats[,1],sim_stats[,6], "SC") ,cbind(sim_stats[,4],sim_stats[,6], "SV") , cbind(sim_stats[,5],sim_stats[,6], "Total") )
#
# #
# colnames(plot_tab)<- c("Recall", "Coverage", "Type")
# plot_tab<- as.data.frame(plot_tab)
# plot_tab$Recall <- as.numeric(plot_tab$Recall)
# plot_tab$rep <- as.factor(rep(c(1:3),15))
# plot_tab_HW<- plot_tab

sim_stats_av
# plot_tab_HW$Coverage<- factor(plot_tab_HW$Coverage, levels = c("3x","5x","8x","12x", "20x"))

sim_stats_av$Coverage<- factor(sim_stats_av$Coverage, levels= c(3,5,8,12,20))

S_sv<-ggplot(sim_stats_av, aes(x = Samplesize, y = M_Correct_SV, color = Coverage, group = Coverage)) +
  geom_point(size = 5, alpha=0.3) +
  geom_line(size=1.1) +
  geom_errorbar(aes(ymin = M_Correct_SV - SD_Correct_SV, ymax = M_Correct_SV +SD_Correct_SV), width = 2, size=1.1) +  # Adjust 'SD' to your column name
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(12.5, 102.5), breaks = c(15, 20, 50, 100)) +
  geom_segment(y = 0, yend = 1, x = 8, xend = 8, color = "black") +
  geom_segment(y = -0.05, yend = -0.05, x = 15, xend = 100, color = "black") +
  labs(x = "", y = "", color = "") +
  scale_color_manual(
    breaks = c(3, 5, 8, 12, 20),
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  theme(
    axis.text.x = element_text(size = 36, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5)
  ) +
  guides(
    shape = "none",
    color = "none"
  )
head(sim_stats_av)
S_sc<-ggplot(sim_stats_av, aes(x = Samplesize, y = M_Correct_SC, color = Coverage, group = Coverage)) +
  geom_point(size = 5, alpha=0.3) +
  geom_line(size=1.1) +
  geom_errorbar(aes(ymin = M_Correct_SC - SD_Correct_SC, ymax = M_Correct_SC +SD_Correct_SC), width = 2, size=1.1) +  # Adjust 'SD' to your column name
  scale_y_continuous(limits = c(0.985, 1)) +
  scale_x_continuous(limits = c(12.5, 102.5), breaks = c(15, 20, 50, 100)) +
  geom_segment(y = 0.985, yend = 1, x = 8, xend = 8, color = "black") +
  geom_segment(y = 0.98425, yend = 0.98425, x = 15, xend = 100, color = "black") +
  labs(x = "Sample size", y = "", color = "") +
  scale_color_manual(
    breaks = c(3, 5, 8, 12, 20),
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
  ) +
  theme(
    axis.text.x = element_text(size = 36, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5)
  ) +
  guides(
    shape = "none",
    color = "none"
  )

library(cowplot)
plot_grid(plotlist = list(HW_sv,S_sv,HW_sc,S_sc ), nrow=2, align = "tblr", rel_widths = c(0.5,0.5))



ggplot(sim_stats_av, aes(x = Samplesize, y = M_Correct_SC, color = Coverage, group = Coverage)) +
  geom_point(size = 8, alpha=1) +
  geom_line(size=3) +
  geom_errorbar(aes(ymin = M_Correct_SC - SD_Correct_SC, ymax = M_Correct_SC +SD_Correct_SC), width = 2, size=1.1) +  # Adjust 'SD' to your column name
  scale_y_continuous(limits = c(0.985, 1)) +
  scale_x_continuous(limits = c(12.5, 102.5), breaks = c(15, 20, 50, 100)) +
  geom_segment(y = 0.985, yend = 1, x = 8, xend = 8, color = "black") +
  geom_segment(y = 0.98425, yend = 0.98425, x = 15, xend = 100, color = "black") +
  labs(x = "Sample size", y = "", color = "Coverage") +
  scale_color_manual(
    breaks = c(3, 5, 8, 12, 20),
    values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
    labels = c("3x", "5x", "8x", "12x", "20x")
  ) +
  theme(
    axis.text.x = element_text(size = 36, angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 40),
    panel.background = element_rect(fill = NA, color = "white"),
    axis.title = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.key.size = unit(2, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 40),
    axis.ticks.length.y = unit(0.5, "cm"),
    axis.ticks.length.x = unit(0.5, "cm"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 40, hjust = 0.5, vjust = 1.5)
  ) +
  guides(
    shape = "none")
