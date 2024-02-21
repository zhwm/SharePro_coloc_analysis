library(ggplot2)

get_H4 <- function(dir,ite=50){
  H4v <- c()
  for (x in 1:ite){
    cs <- read.table(paste0(dir,'SH10_1e-3/','Q',x,'.z_G',x,'.z.cs'),header=T)
    H4 <- ifelse(dim(cs)[1],max(cs$share),0)
    H4v <- c(H4v, H4)
  }
  return(H4v)
}

get_coloc_H4 <- function(dir){
  coloc <- read.table(paste0(dir,'coloc_1e-3.txt'),header=T)
  return(coloc$PP.H4.abf)
}

wdir = '../sim/'

SHv <- c()
colocv <- c()
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
      colocv <- c(colocv,get_coloc_H4(dir))
      lov <- c(lov,rep(lo,50))
      kcv <- c(kcv,rep(kc,50))
      ksv <- c(ksv,rep(ks,50))
    }
  }
}

dat = data.frame(SharePro=SHv, COLOC=colocv, KC=kcv, KS=ksv, LOCI=lov)
write.table(dat,"../doc/sharepro_loc_sim_colocalization_1e-3.csv",sep='\t', row.names = F, quote = F)

