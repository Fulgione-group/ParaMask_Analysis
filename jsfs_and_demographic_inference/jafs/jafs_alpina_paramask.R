
# Run as:
# Rscript jsfs_cvi_toPublish.R

###
#
# Set working directory to the CVI directory from github
#
###
setwd("/Users/fulgione/Documents/useful_stuff/alpina/jsfs/")





# Load required libraries 
library("spam")
library("fields")

par(ps=7)
colors <- c("white", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000")






jsfs = as.matrix(read.table("Path to jsfs file produced in java_programs", head=F))

jsfs[nrow(jsfs), ncol(jsfs)] = 0
jsfs[1,1] = 0

# All SNPs:
allSnps = sum(jsfs)
allSnps


# Segregating in Pop1:
segS=0
for (row in 2:(nrow(jsfs)-1)) {
  for (col in 1:ncol(jsfs)) {
    segS = segS + jsfs[row,col]
  }
}
segS

# Segregating in Pop2:
segF=0
for (col in 2:(ncol(jsfs)-1)) {
  for (row in 1:nrow(jsfs)) {
    segF = segF + jsfs[row,col]
  }
}
segF

# Shared polymorphisms
sharedPolim=0
for (col in 2:(ncol(jsfs)-1)) {
  for (row in 2:(nrow(jsfs)-1)) {
    sharedPolim = sharedPolim + jsfs[row,col]
  }
}
sharedPolim

# % shared variation:
sharedPolim/allSnps

jsfs = jsfs + 1
jsfs=log10(jsfs)

ticks<- c(0, 10, 100, 1000, 10000, 100000, 1000000)
labelS=c(expression(0), expression(10), expression(10^"2"), expression(10^"3"), expression(10^"4"), expression(10^"5"), expression(10^"6"))

pdf("Path to pdf", width=2.2, height=1.6)
par(ps=7, mfrow=c(1,1), mar=c(1.1,1.1,0.5,0.5))

image.plot(t(jsfs), axis.args=list(line=-0.8, at=c(0, log10(ticks[2:length(ticks)])), labels=labelS, lwd=-0, lwd.ticks=-0.5, tck=-0.5), col=colors, axes=F, legend.width=0.5, smallplot = c(.8, .84, .2, .8))

title(main="")
axis(2, at=c(-01, 1.1), tick=T, labels=c("a", "a"), lwd.ticks=-1, lwd=1)          # , cex.axis=1.5
axis(2, las=2, at=c(0, 1), tick=F, labels=c("0", dim(jsfs)[[1]]-1), line=-0.9)
axis(1, at=c(-01, 1.1), tick=T, labels=c("a", "a"), lwd.ticks=-1, lwd=1)
axis(1, at=c(0, 1), tick=F, labels=c("0", dim(jsfs)[[2]]-1), line=-1.3)
axis(4, at=c(-01, 1.1), tick=T, labels=c("a", "a"), lwd.ticks=-1, lwd=1)
axis(3, at=c(-01, 1.1), tick=T, labels=c("a", "a"), lwd.ticks=-1, lwd=1)

title("", xlab=expression(paste("ES04")), ylab=expression(paste("ES03")), line=0.2)

dev.off()

















###
###
###
#
#     Let us look at the SFS
#
###
###
###

#
sfs1 = as.matrix(read.table("Path to sfs1", head=F))
sfs2 = as.matrix(read.table("Path to sfs2", head=F))

sfs1 = sfs1[2:(length(sfs1)-1)] / sum(sfs1[2:(length(sfs1)-1)])
sfs2 = sfs2[2:(length(sfs2)-1)] / sum(sfs2[2:(length(sfs2)-1)])
twosfss<-rbind(sfs1, sfs2)

pdf("Path to pdf", width=2.2, height=1.6)
par(ps=7, mfrow=c(1,1), mar=c(1.1,1.1,0.5,0.5))

mp <- barplot(twosfss, main="", space=c(0.1, 0.4),
        xlab="", col=c("goldenrod1", "darkgoldenrod4"), # darkblue","red"),
        #legend = rownames(counts), 
        beside=TRUE, axes=F, ylim=c(0, 0.2), border=NA, names.arg=rep("", 29))

mtext(side = 1, at = c(colMeans(mp)[[1]], colMeans(mp)[[29]]), line = -0.3,
      text = c(1, 29), col = "black")

title(main="")
axis(1, at=c(colMeans(mp)-1.1, colMeans(mp)[[length(colMeans(mp))]] + 1.1), tick=T, labels=F, lwd=0.4, line=0.04, tck=-0.01)
# axis(1, at=c(colMeans(mp)-1.1, colMeans(mp)[[length(colMeans(mp))]] + 1.1), tick=F, labels=T, lwd=0.4, line=0.04, tck=-0.01)
axis(2, las=2, at=seq(0, 0.2, 0.05), tick=T, labels=F, line=-0.2, lwd=0.4, tck=-0.01)
axis(2, las=2, at=seq(0, 0.2, 0.05), tick=F, labels=seq(0, 20, 5), line=-1.1, lwd=0.4, tck=-0.01)

title("", xlab=expression(paste("Number of derived alleles (n=30)")), line=0.2)
title("", ylab=expression(paste("Frequency "(10^"-2"))), line=0.3)

legend("topright", 
       legend = c("ES03", "ES04"), 
       fill = c("goldenrod1", "darkgoldenrod4"), border=NA, bty="n")

dev.off()

