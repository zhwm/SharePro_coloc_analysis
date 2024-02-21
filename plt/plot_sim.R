library(ggplot2)
setwd('~/Desktop/Projects/SharePro_loc/plt')

dat = read.table("../doc/sharepro_loc_sim_colocalization.csv",sep='\t',header=T)

pltdat = data.frame(Prob=c(dat$SharePro, dat$CSuSiE, dat$PWCoCo, dat$COLOC, dat$eCAVIAR),
                    method=rep(c("SharePro", "COLOC+SuSiE", "PWCoCo", "COLOC", "eCAVIAR"),each=nrow(dat)),
                     KC=rep(dat$KC,5),
                     KS=rep(dat$KS,5)) 
pltdat$H4 = ifelse(pltdat$KC!=0,"Colocalized","Non-colocalized")
pltdat$H4 = factor(pltdat$H4, levels = rev(c("Non-colocalized","Colocalized")))
pltdat$K = paste0("K[C]~+~K[S]:",pltdat$KC + pltdat$KS)

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
ggsave("../doc/SharePro_loc_sim.pdf",height = 15, width = 15)


dat8 = aggregate(Prob>0.8 ~ KC + KS + method, dat=pltdat, mean)
dat8summary = reshape(dat8, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat8summary$'KC+KS' = dat8summary$KC + dat8summary$KS
names(dat8summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(dat8summary,"../doc/sharepro_loc_sim_0.8.csv",sep=',', row.names = F, quote = F)

dat6 = aggregate(Prob>0.6 ~ KC + KS + method, dat=pltdat, mean)
dat6summary = reshape(dat6, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat6summary$'KC+KS' = dat6summary$KC + dat6summary$KS
names(dat6summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(dat6summary,"../doc/sharepro_loc_sim_0.6.csv",sep=',', row.names = F, quote = F)

dat4 = aggregate(Prob>0.4 ~ KC + KS + method, dat=pltdat, mean)
dat4summary = reshape(dat4, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat4summary$'KC+KS' = dat4summary$KC + dat4summary$KS
names(dat4summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(dat4summary,"../doc/sharepro_loc_sim_0.4.csv",sep=',', row.names = F, quote = F)

dat2 = aggregate(Prob>0.2 ~ KC + KS + method, dat=pltdat, mean)
dat2summary = reshape(dat2, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat2summary$'KC+KS' = dat2summary$KC + dat2summary$KS
names(dat2summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro","KC+KS")
write.table(dat2summary,"../doc/sharepro_loc_sim_0.2.csv",sep=',', row.names = F, quote = F)

datn8 = aggregate(Prob<0.8 ~ KC + KS + method, dat=pltdat, mean)
datn8summary = reshape(datn8, direction = "wide", timevar = "method", idvar = c('KS','KC'))
datn8summary$'KC+KS' = datn8summary$KC + datn8summary$KS
names(datn8summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(datn8summary,"../doc/sharepro_loc_sim_0.8_neg.csv",sep=',', row.names = F, quote = F)

datn6 = aggregate(Prob<0.6 ~ KC + KS + method, dat=pltdat, mean)
datn6summary = reshape(datn6, direction = "wide", timevar = "method", idvar = c('KS','KC'))
datn6summary$'KC+KS' = datn6summary$KC + datn6summary$KS
names(datn6summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(datn6summary,"../doc/sharepro_loc_sim_0.6_neg.csv",sep=',', row.names = F, quote = F)

datn4 = aggregate(Prob<0.4 ~ KC + KS + method, dat=pltdat, mean)
datn4summary = reshape(datn4, direction = "wide", timevar = "method", idvar = c('KS','KC'))
datn4summary$'KC+KS' = datn4summary$KC + datn4summary$KS
names(datn4summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(datn4summary,"../doc/sharepro_loc_sim_0.4_neg.csv",sep=',', row.names = F, quote = F)

datn2 = aggregate(Prob<0.2 ~ KC + KS + method, dat=pltdat, mean)
datn2summary = reshape(datn2, direction = "wide", timevar = "method", idvar = c('KS','KC'))
datn2summary$'KC+KS' = datn2summary$KC + datn2summary$KS
names(datn2summary) = c("KC","KS","COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro", "KC+KS")
write.table(datn2summary,"../doc/sharepro_loc_sim_0.2_neg.csv",sep=',', row.names = F, quote = F)


## time table
timedt = read.csv('../doc/time.summary',header=F,sep=' ')
timedt$method=c("eCAVIAR",
                "COLOC+SuSiE",
                "COLOC",
                "PWCoCo",
                "SharePro")
dftime = aggregate(V2~method,data=timedt,mean)
dftime$sd = aggregate(V2~method,data=timedt,sd)$V2
colnames(dftime) = c("Method","Run time average(s)", "Run time standard deviation(s)")
dftime$`Run time average(s)` = dftime$`Run time average(s)`/50
dftime$`Run time standard deviation(s)`= dftime$`Run time standard deviation(s)`/50
write.table(dftime,'../doc/time.csv',sep=',',col.names = T,row.names = F)

## Sensitivity analyses
sdir = '../sim/example/'
prior = c('1e-3','5e-4','2e-4','1e-4','5e-5','2e-5','1e-5','5e-6','2e-6','1e-6','5e-7','2e-7','1e-7')
ex = c(14,15)
colprob = c()

for (pr in prior) {
  for (x in ex) {
    cs = read.table(paste0(sdir,pr,'/Q',x,'.fastGWA.z_G',x,'.fastGWA.z.cs'),header=T)
    colprob = c(colprob,max(cs$share))
  }
}

sdf = data.frame(coloc = colprob,
                  prior = rep(prior,each=2),
                  ex = rep(ex,length(prior)))
sdf$prior = as.numeric(sdf$prior)
ggplot(sdf[sdf$ex==15,],aes(x=log10(prior),y=coloc)) + geom_point()  + theme_bw() + ylab("Colocalization probability") +
  scale_y_continuous(breaks = c(0,0.2,0.8,1), limits = c(0,1)) +
  xlab(expression(log[10](prior))) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13))
ggsave("../doc/prior_sensitivity.pdf",height = 3, width = 5)

ggplot(sdf[sdf$ex==14,],aes(x=log10(prior),y=coloc)) + geom_point()  + theme_bw() +   ylab("Colocalization probability") + 
  scale_y_continuous(breaks = c(0,0.2,0.8,1), limits = c(0,1)) +
  xlab(expression(log[10](prior))) + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13))
ggsave("../doc/prior_insensitivity.pdf",height = 3, width = 5)

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
  ggsave(paste0("../doc/prior_sensitivity",x,".pdf"),height = 5, width = 5)
}

gwasplot(14,"rs175860")
gwasplot(15,"rs12205824")
