library(ggplot2)

get_H4 <- function(dir,ite=50){
  H4v <- c()
  for (x in 1:ite){
    cs <- read.table(paste0(dir,'SH/','Q',x,'.z_G',x,'.z.cs'),header=T)
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

wdir = '~/scratch/SharePro_loc/sim/'

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
write.table(dat,"~/scratch/SharePro_loc/doc/sharepro_loc_sim_colocalization.csv",sep='\t', row.names = F, quote = F)

dat <- read.table("~/scratch/SharePro_loc/doc/sharepro_loc_sim_colocalization.csv",sep='\t',header=T)

pltdat <- data.frame(Prob=c(dat$SharePro,
                            dat$CSuSiE,
                            dat$PWCoCo,
                            dat$COLOC,
                            dat$eCAVIAR),
                    method=rep(c("SharePro",
                                 "COLOC+SuSiE",
                                 "PWCoCo",
                                 "COLOC",
                                 "eCAVIAR"),each=nrow(dat)),
                    KC=rep(dat$KC,5),
                    KS=rep(dat$KS,5)) 
pltdat$H4 <- ifelse(pltdat$KC!=0,"Colocalized","Non-colocalized")
pltdat$H4 <- factor(pltdat$H4, levels = rev(c("Non-colocalized","Colocalized")))
pltdat$K <- paste0("K[C]~+~K[S]:",pltdat$KC + pltdat$KS)

ggplot(dat=pltdat,aes(x=method,y=Prob,color = as.factor(KC))) + geom_boxplot(position = position_dodge(width = 0.6), width = 0.3) + 
  theme_bw() + xlab("") + 
  theme(axis.title = element_text(size = 19),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14),
        legend.position = "top") + facet_grid(K~H4,labeller = label_parsed) +
  labs(color = "Number of shared causal variants") +
  scale_color_manual(values = c("orange","brown","red","pink","royalblue","darkblue")) +
  ylab("Colocalization probability") +
  scale_y_continuous(breaks = c(0.2,0.8))
ggsave("~/scratch/SharePro_loc/doc/sharepro_loc_sim_colocalization.pdf",height = 15, width = 15)

# summary table
dat8 = aggregate(Prob>0.8 ~ KC + KS + method, dat=pltdat, mean)
dat8summary = reshape(dat8, direction = "wide", timevar = "method", idvar = c('KS','KC'))
names(dat8summary) <- c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro")
write.table(dat8summary,"~/scratch/SharePro_loc/doc/sharepro_loc_sim_0.8.csv",sep=',', row.names = F, quote = F)
dat2 = aggregate(Prob<0.2 ~ KC + KS + method, dat=pltdat, mean)
dat2summary = reshape(dat2, direction = "wide", timevar = "method", idvar = c('KS','KC'))
names(dat2summary) <- c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro")
write.table(dat2summary,"~/scratch/SharePro_loc/doc/sharepro_loc_sim_0.2.csv",sep=',', row.names = F, quote = F)

## time table
timedt = read.csv('~/scratch/SharePro_loc/doc/time.summary',header=F,sep=' ')
timedt$method=c("eCAVIAR",
             "COLOC+SuSiE",
             "COLOC",
             "PWCoCo",
             "SharePro")
dftime = aggregate(V2~method,data=timedt,mean)
dftime$sd = aggregate(V2~method,data=timedt,sd)$V2
colnames(dftime) <- c("Method","Run time average(s)", "Run time standard deviation(s)")
dftime$`Run time average(s)` = dftime$`Run time average(s)`/50
dftime$`Run time standard deviation(s)`= dftime$`Run time standard deviation(s)`/50
write.table(dftime,'~/scratch/SharePro_loc/doc/time.csv',sep=',',col.names = T,row.names = F)

## Sensitivity analyses
sdir = '~/scratch/SharePro_loc/sim/example/'
prior <- c('1e-3','5e-4','2e-4','1e-4','5e-5','2e-5','1e-5','5e-6','2e-6','1e-6','5e-7','2e-7','1e-7')
ex <- c(14,15)
colprob <- c()

for (pr in prior) {
  for (x in ex) {
    cs <- read.table(paste0(sdir,pr,'/Q',x,'.fastGWA.z_G',x,'.fastGWA.z.cs'),header=T)
    colprob <- c(colprob,max(cs$share))
  }
}

sdf <- data.frame(coloc = colprob,
                  prior = rep(prior,each=2),
                  ex = rep(ex,length(prior)))
sdf$prior <- as.numeric(sdf$prior)
ggplot(sdf[sdf$ex==15,],aes(x=log10(prior),y=coloc)) + geom_point()  + theme_bw() + ylab("Colocalization probability") +
  scale_y_continuous(breaks = c(0,0.2,0.8,1), limits = c(0,1)) +
  xlab(expression(log[10](prior))) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13))
ggsave("~/scratch/SharePro_loc/doc/prior_sensitivity.pdf",height = 3, width = 5)
ggplot(sdf[sdf$ex==14,],aes(x=log10(prior),y=coloc)) + geom_point()  + theme_bw() +   ylab("Colocalization probability") + 
  scale_y_continuous(breaks = c(0,0.2,0.8,1), limits = c(0,1)) +
  xlab(expression(log[10](prior))) + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13))
ggsave("~/scratch/SharePro_loc/doc/prior_insensitivity.pdf",height = 3, width = 5)

gwasplot <- function(x,causal){
  gwas <- read.table(paste0(sdir,'G',x,'.fastGWA'),header=T)
  qtl <- read.table(paste0(sdir,'Q',x,'.fastGWA'),header=T)
  ld <- read.table(paste0(sdir,'Locus4.ld'))
  pltdat <- data.frame(Pos = rep(gwas$POS,2),
                       p = c(gwas$P,qtl$P),
                       r2 = rep((ld[,which(gwas$SNP==causal)])^2,2),
                       study = rep(c("GWAS of trait 1","GWAS of trait 2"),each=nrow(gwas)))
  ggplot(pltdat, aes(x = Pos / 1000000, y = -log10(p), color = r2)) +
    geom_point() +
    facet_grid(study~., scale = "free") +
    scale_color_stepsn(breaks = c(0,0.2,0.4,0.6,0.8,1), colors = c("darkblue","royalblue","green","orange","red")) +
    theme_bw() + 
    ylab(expression(-log[10](p-value))) +
    xlab("Chr6 (Mb)") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13)) +
    labs(color = expression(r^2))
  ggsave(paste0("~/scratch/SharePro_loc/doc/prior_sensitivity",x,".pdf"),height = 5, width = 5)
}

gwasplot(14,"rs175860")
gwasplot(15,"rs12205824")
