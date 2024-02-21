library(coloc)

args = commandArgs(trailingOnly=T)
gf=args[1]
qf=args[2]
lf=args[3]
prior=args[4]
out=args[5]

dat <- data.frame(nsnps = numeric(),hit1 = character(),hit2 = character(), PP.H0.abf = numeric(),PP.H1.abf=numeric(),PP.H2.abf=numeric(),PP.H3.abf=numeric(),PP.H4.abf=numeric(),idx1=numeric(),idx2=numeric())

gwas <- read.table(gf, header=T)
qtl <- read.table(qf, header=T)
LD <- as.matrix(read.table(lf))
rownames(LD) <- colnames(LD) <- qtl$SNP
q <- list(snp=qtl$SNP,beta=qtl$BETA,varbeta=qtl$SE**2,type="quant",MAF=as.numeric(qtl$EAF),N=qtl$N,LD=LD)
g <- list(snp=gwas$SNP,beta=gwas$BETA,varbeta=gwas$SE**2,type="quant",MAF=as.numeric(gwas$EAF),N=gwas$N,LD=LD)
sqtl <- runsusie(q)
sgwas <- runsusie(g)
res <- coloc.susie(q,g,p12=as.numeric(prior))
output <- as.data.frame(res$summary)
print(dim(output)[1])
if (dim(output)[1]>0){
  dat <- rbind.data.frame(dat,output)
}
write.table(dat, out, col.names=T,row.names=F,sep="\t",quote=F)
