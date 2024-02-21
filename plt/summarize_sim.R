library(ggplot2)

get_H4 <- function(dir,ite=50){
  H4v <- c()
  for (x in 1:ite){
    cs <- read.table(paste0(dir,'SH10/','Q',x,'.z_G',x,'.z.cs'),header=T)
    H4 <- ifelse(dim(cs)[1],max(cs$share),0)
    H4v <- c(H4v, H4)
  }
  return(H4v)
}

get_csusie_H4 <- function(dir,ite=50){
  H4v <- c()
  for (x in 1:ite){
    H4 <- 0
    resf <- paste0(dir,x,'_coloc_susie.txt')
    if (file.exists(resf)){
      cs <- read.table(resf,header=T)$PP.H4.abf
      H4 <- ifelse(length(cs)>0,max(cs),0)
    }
    H4v <- c(H4v, H4)
  }
  return(H4v)
}

get_pwcoco_H4 <- function(dir){
  pwcoco = read.table(paste0(dir,'pwcoco_out.coloc'),header=T)
  aggregate(H4 ~ Dataset1, pwcoco, max)$H4
}

get_coloc_H4 <- function(dir){
  coloc <- read.table(paste0(dir,'coloc.txt'),header=T)
  return(coloc$PP.H4.abf)
}

get_ecaviar_H4 <- function(dir,ite=50){
  H4v <- c()
  for (x in 1:ite){
    ecaviar <- read.table(paste0(dir,x,'_col'),header=T)
    H4 <- max(ecaviar$CLPP)
    H4v <- c(H4v, H4)
  }
  return(H4v)
}

wdir = '../sim/'

SHv <- c()
csusiev <- c()
pwcocov <- c()
colocv <- c()
ecav <- c()
lov <- c()
kcv <- c()
ksv <- c()

for (l in 1:5) {
  for (k in 1:5) {
    for (kc in 0:k) {
      lo = paste0('Locus',l,'/')
      ks = k - kc
      dir = paste0(wdir,lo,paste0(kc,"_", ks,"_0.05_0.01/"))
      print(dir)
      SHv <- c(SHv,get_H4(dir))
      csusiev <- c(csusiev,get_csusie_H4(dir))
      pwcocov <- c(pwcocov,get_pwcoco_H4(dir))
      colocv <- c(colocv,get_coloc_H4(dir))
      ecav <- c(ecav,get_ecaviar_H4(dir))
      lov <- c(lov,rep(lo,50))
      kcv <- c(kcv,rep(kc,50))
      ksv <- c(ksv,rep(ks,50))
    }
  }
}

dat = data.frame(SharePro=SHv,
                 CSuSiE=csusiev,
                 PWCoCo=pwcocov,
                 COLOC=colocv,
                 eCAVIAR = ecav,
                 KC=kcv,
                 KS=ksv,
                 LOCI=lov)
write.table(dat,"../doc/sharepro_loc_sim_colocalization.csv",sep='\t', row.names = F, quote = F)

