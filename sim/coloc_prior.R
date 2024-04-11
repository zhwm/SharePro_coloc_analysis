library(coloc)
args = commandArgs(trailingOnly=T)
prior=args[1]

dat <- matrix(nrow = 50, ncol = 6)

for (ite in 1:50){
  qtl <- read.table(paste0("Q",ite,".fastGWA"),header = T)
  gwas <- read.table(paste0("G",ite,".fastGWA"),header = T)
  q <- list(snp=qtl$SNP,beta=qtl$BETA,varbeta=qtl$SE**2,position=qtl$POS,type="quant",MAF=qtl$AF1,N=qtl$N)
  g <- list(snp=gwas$SNP,beta=gwas$BETA,varbeta=gwas$SE**2,position=gwas$POS,type="quant",MAF=gwas$AF1,N=gwas$N)
  abf.res <- coloc.abf(dataset1=q,dataset2=g,p12=as.numeric(prior))
  dat[ite,] <- as.vector(abf.res$summary)
  print(abf.res$summary)
}
colnames(dat) <- names(abf.res$summary)
write.table(dat, paste0("coloc_",prior,".txt"), col.names=T,row.names=F,sep="\t",quote=F)
