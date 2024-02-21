library(ggplot2)
setwd('~/Desktop/Projects/SharePro_loc/plt')

dat = read.table("../doc/sharepro_loc_sim_colocalization_1e-3.csv",sep='\t',header=T)

pltdat = data.frame(Prob=c(dat$SharePro, dat$COLOC),
                    method=rep(c("SharePro", "COLOC"),each=nrow(dat)),
                    KC=rep(dat$KC,2),
                    KS=rep(dat$KS,2)) 
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
ggsave("../doc/SharePro_loc_sim_1e-3.pdf",height = 15, width = 15)


dat8 = aggregate(Prob>0.8 ~ KC + KS + method, dat=pltdat, mean)
dat8summary = reshape(dat8, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat8summary$'KC+KS' = dat8summary$KC + dat8summary$KS
names(dat8summary) = c("KC","KS","COLOC","SharePro", "KC+KS")
write.table(dat8summary,"../doc/sharepro_loc_sim_0.8_1e-3.csv",sep=',', row.names = F, quote = F)

dat2 = aggregate(Prob<0.2 ~ KC + KS + method, dat=pltdat, mean)
dat2summary = reshape(dat2, direction = "wide", timevar = "method", idvar = c('KS','KC'))
dat2summary$'KC+KS' = dat2summary$KC + dat2summary$KS
names(dat2summary) = c("KC","KS","COLOC","SharePro", "KC+KS")
write.table(dat2summary,"../doc/sharepro_loc_sim_0.2_1e-3.csv",sep=',', row.names = F, quote = F)

pwer = read.table("../doc/sharepro_loc_sim_cs_power.csv",sep=',',header=T)
plotdat = data.frame(Prob=c(pwer$COLOC.SuSiE, pwer$eCAVIAR, pwer$SharePro),
                     method = rep(c("COLOC+SuSiE", "eCAVIAR", "SharePro"), each=nrow(pwer)),
                     KC = rep(pwer$KC,3), 
                     KS = rep(pwer$KS, 3))
pwermean = aggregate(Prob ~ KC + KS + method, dat=plotdat, mean)
pwersd = aggregate(Prob ~ KC + KS + method, dat=plotdat, sd)
pwermean$sd = pwersd$Prob
pwermean$K = paste0("K[C]~+~K[S]:",pwermean$KC + pwermean$KS)
pwermean$upper = pwermean$Prob + 1.96 * pwermean$sd
pwermean$upper[pwermean$upper > 1] = 1
pwermean$lower = pwermean$Prob - 1.96 * pwermean$sd
pwermean$lower[pwermean$lower < 0] = 0
ggplot(dat=pwermean, aes(x=method, y=Prob, color=as.factor(KC))) + geom_point(position = position_dodge(width = 0.6)) + facet_grid(K~., labeller = label_parsed) + 
  scale_color_manual(values = c("brown","red","pink","royalblue","darkblue")) + ylab("Power") +
  geom_errorbar(aes(ymin=lower, ymax=upper), position = position_dodge(width = 0.6), width=0.1) + theme_bw() + xlab("") +
  labs(color = "Number of shared causal variants") + 
  theme(axis.title = element_text(size = 19),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = "top")
ggsave("../doc/SharePro_loc_cs_power.pdf",height = 15, width = 15)


