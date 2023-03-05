args = commandArgs(trailingOnly=TRUE) #LOCI KC KS hq hg
set.seed(0)
library(data.table)

LO = args[1]
KC = as.numeric(args[2])
KS = as.numeric(args[3])
varq = as.numeric(args[4])
varg = as.numeric(args[5])

fq = fread(paste0(LO,"_Q.raw"))
fg = fread(paste0(LO,"_G.raw"))

for (j in 7:ncol(fq)) {
  if (colnames(fq)[j]!=colnames(fg)[j]) {
    fq[,7] = 2 - fq[,7]
    colnames(fq)[j] = colnames(fg)[j]
  }
}

allyq = list(fq$FID,fq$IID)
allyg = list(fg$FID,fg$IID)
allcsnp = list()
ite=1
while (length(allyq)<52) {
  csnp = colnames(fq)[sample(7:ncol(fq),(KC+2*KS))]
  genoq = scale(fq[,csnp[1:(KC+KS)],with=F])
  genog = scale(fg[,csnp[(KS+1):(KC+KS*2)],with=F])
  yqhat = as.matrix(genoq) %*% rep(sqrt(varq),(KC+KS))
  yghat = as.matrix(genog) %*% rep(sqrt(varg),(KC+KS))
  yq = yqhat + rnorm(length(yqhat),0,sqrt(1-var(yqhat,na.rm = T)))
  yg = yghat + rnorm(length(yghat),0,sqrt(1-var(yghat,na.rm = T)))
  if (max(min(coef(summary(lm(yg~genog)))[-1,4]),min(coef(summary(lm(yq~genoq)))[-1,4]))<1e-5){
    allyq[[ite+2]]=as.vector(yq)
    allyg[[ite+2]]=as.vector(yg)
    allyq[[ite+2]][is.na(allyq[[ite+2]])] <- 'NA'
    allyg[[ite+2]][is.na(allyg[[ite+2]])] <- 'NA'
    allcsnp[[ite]]=unlist(lapply(strsplit(csnp,'_'),function(x){x[1]}))
    ite=ite+1    
  }
}

setDT(allyq)
setDT(allyg)
setnames(allyq, c(c('FID','IID'),paste0('trait',seq(1,ite-1))))
setnames(allyg, c(c('FID','IID'),paste0('trait',seq(1,ite-1))))
fwrite(allyq,file='Q.phen',sep='\t',quote=F,col.names=F) #Pheno file
fwrite(allyg,file='G.phen',sep='\t',quote=F,col.names=F) #Pheno file
write.table(do.call(rbind,allcsnp),file='csnp.txt',row.names = F,col.names=F)
