#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

vcf <- args[1]



library(rCNV)

vcf <- readVCF(vcf, verbose = FALSE)
ad.tab <- hetTgen(vcf, info.type= "AD")
gt.tab <- hetTgen(vcf, info.type = "GT")
ad.tab <- ad.correct(ad.tab, gt.table = gt.tab)

pdf("Rplot.norm.allele.depth.pdf")
ad.nor <- cpm.normal(ad.tab, method="MedR", plot = TRUE)
dev.off()

pdf("Rplot.allele.cov.pdf")
A.info.tab <- allele.info(X = ad.tab, x.norm = ad.nor, plot.allele.cov = TRUE)
dev.off()

pdf("Rplot_deviants.detected.pdf")
deviants <- dupGet(A.info.tab, test = c("z.05","chi.05"), plot = TRUE, verbose = TRUE)
dev.off()

pdf("Rplot_CNVs.pdf")
CV <- cnv(A.info.tab, test=c("z.05","chi.05"), filter = "kmeans", plot = TRUE)
dev.off()


write.table(deviants, file = "deviants.txt")
write.table(CV, file = "cnv.txt")
