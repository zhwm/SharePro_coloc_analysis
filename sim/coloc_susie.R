library(coloc)
args = commandArgs(trailingOnly=T)
region=args[1]
ite=args[2]

dat <- data.frame(nsnps = numeric(),hit1 = character(),hit2 = character(), PP.H0.abf = numeric(),PP.H1.abf=numeric(),PP.H2.abf=numeric(),PP.H3.abf=numeric(),PP.H4.abf=numeric(),idx1=numeric(),idx2=numeric(),ite=numeric())

print(ite)
qtl <- read.table(paste0("Q",ite,".fastGWA"),header = T)
gwas <- read.table(paste0("G",ite,".fastGWA"),header = T)
qtlLD <- as.matrix(read.table(paste0(region,".ld")))
rownames(qtlLD) <- colnames(qtlLD) <- qtl$SNP
q <- list(snp=qtl$SNP,beta=qtl$BETA,varbeta=qtl$SE**2,position=qtl$POS,type="quant",MAF=qtl$AF1,N=qtl$N,LD=qtlLD)
g <- list(snp=gwas$SNP,beta=gwas$BETA,varbeta=gwas$SE**2,position=gwas$POS,type="quant",MAF=gwas$AF1,N=gwas$N,LD=qtlLD)
sqtl <- runsusie(q)
sgwas <- runsusie(g)
res <- coloc.susie(q,g)
output <- as.data.frame(res$summary)
print(dim(output)[1])
if (dim(output)[1]>0){
  output$ite <- ite
  dat <- rbind.data.frame(dat,output)
}
write.table(dat, paste0(ite,"_coloc_susie.txt"), col.names=T,row.names=F,sep="\t",quote=F)
