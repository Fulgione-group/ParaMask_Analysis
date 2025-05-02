###read siumulation mat file

library(MASS)

args = commandArgs(trailingOnly=TRUE)

matpath<- args[1]
varcov <- as.numeric(args[2])
outpath <- paste(paste(strsplit(matpath,split = "\\.")[[1]][1:(length(strsplit(matpath,split = "\\.")[[1]])-1)] , collapse = "."),"_cov", args[2] , ".vcf", sep="")

mat<- read.table(file = matpath, header = F, sep = "\t")



# head(mat)


  

header1 <- "##fileformat=VCFv4.2"
header2<-  c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste0("genotype_", 1:100))
header2 <-as.data.frame(t(header2))

infofields<- c("chr1", "pos", "ID", "T", "A", ".", ".", "DP=372;AN=4", "GT:AD:DP:RGQ")



convert_snp<- function(snp, rho){
  geno<-c()
  SVinfo<-as.character(snp[length(snp)])
  tmp_info<- infofields
  tmp_info[2] <- as.numeric(snp[1])
  tmp_info[3] <-paste(tmp_info[1],tmp_info[2],SVinfo, sep = "_")[[1]][1]
  if(startsWith(as.character(snp[length(snp)]),"SV")){
    SV <- TRUE
    cov <- rnegbin(n = (length(snp)-3), mu = varcov*2, theta =8)
  }else{
    SV<-FALSE
    cov <- rnegbin(n = (length(snp)-3), mu = varcov, theta = 8)
  }
  cov[cov<0] <- 0
  for(i in 2:(length(snp)-2)){
    if(cov[(i-1)]<=1){
      geno<- c(geno, paste("./.",".,.", cov[(i-1)],"40", sep=":"))
    }else{
      if(snp[i]==0){
        geno<- c(geno, paste("0/0",paste(cov[(i-1)], "0",sep=","), cov[(i-1)],"40", sep=":"))
      }else if(snp[i]==2){
        geno<- c(geno, paste("1/1",paste("0",cov[(i-1)],sep=","), cov[(i-1)],"40", sep=":"))
      } else if(snp[i]==1){
        # ad<-rbetabinom(n = 1,size=cov[(i-1)], prob =0.5 , rho = rho )
        ad<-rbinom(n = 1, size=cov[(i-1)], prob =0.5)
        ar<-ad/cov[(i-1)]
        if(ar>0.8){
          geno<- c(geno, paste("1/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else if(ar<0.2){
          geno<- c(geno, paste("0/0",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else{
          geno<- c(geno, paste("0/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }
      } else if(snp[i]==3){
        # ad<-rbetabinom(n = 1,size=cov[(i-1)], prob =0.25 , rho = rho )
        ad<-rbinom(n = 1, size=cov[(i-1)], prob =0.25)
        ar<-ad/cov[(i-1)]
        if(ar>0.8){
          geno<- c(geno, paste("1/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else if(ar<0.2){
          geno<- c(geno, paste("0/0",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else{
          geno<- c(geno, paste("0/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }
      } else if(snp[i]==4){
        # ad<-rbetabinom(n = 1,size=cov[(i-1)], prob =0.75 , rho = rho )
        ad<-rbinom(n = 1, size=cov[(i-1)], prob =0.75)
        ar<-ad/cov[(i-1)]
        if(ar>0.8){
          geno<- c(geno, paste("1/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else if(ar<0.2){
          geno<- c(geno, paste("0/0",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }else{
          geno<- c(geno, paste("0/1",paste((cov[(i-1)]-ad),ad,sep=","), cov[(i-1)],"40", sep=":"))
        }
      }
    }  
  }
  line <- as.data.frame(c(tmp_info,geno))
  return(line)
}

vcf_list<- apply(mat, MARGIN = 1, FUN = function(x){convert_snp(snp = x,rho = rho)})

vcf<- as.data.frame(do.call(cbind, vcf_list))
vcf<-as.data.frame(t(vcf))

cat(header1, "\n", file = outpath)
write.table(x = header2, file = outpath, sep="\t", quote = F, row.names = F, col.names = F, append = T)
write.table(x = vcf, file = outpath, sep="\t", quote = F, row.names = F, col.names = F, append = T)




