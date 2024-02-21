library(ggplot2)
setwd('~/Desktop/Projects/SharePro_loc/plt')

gwas = read.table('../dat/RSPO3_eBMD/BMD_pwcoco.txt', header=T)
qtl = read.table('../dat/RSPO3_eBMD/RSPO3_pwcoco.txt', header=T)
cs = read.table('../dat/RSPO3_eBMD/SH_BMD_RSPO3_1e-5', header=T)

gwas$mlogP = -log10(exp(1)) * pchisq((gwas$BETA/gwas$SE)^2, df=1, lower.tail = F, log.p=T)
qtl$mlogP = -log10(exp(1)) * pchisq((qtl$BETA/qtl$SE)^2, df=1, lower.tail = F, log.p=T)

plotdat = data.frame(SNP = gwas$SNP, gwas = gwas$mlogP, qtl = qtl$mlogP, cs = "others")
for (i in 1:nrow(cs)) {
  cslist = strsplit(cs$cs[i], '/')[[1]]
  plotdat$cs[plotdat$SNP %in% cslist] = paste("Effect group", i)
}

plotdat = plotdat[order(plotdat$cs, decreasing = T),]

ggplot(plotdat, aes(x = qtl, y = gwas, color = cs)) +
  geom_vline(xintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size=2) + theme_bw() +
  scale_color_manual(values = c("hotpink", "darkmagenta", "cyan", "cyan4", "coral", "coral4", "grey")) +
  xlab(expression (-log[10](p-value)~"in"~RSPO3~pQTL)) +
  ylab(expression (-log[10](p-value)~"in"~BMD~GWAS)) +
  labs(color = '') +
  scale_size_manual(values = c(5, 2)) +
  theme(axis.title = element_text (size = 14),
        axis.text = element_text (size = 13),
        legend.text = element_text (size = 13))
#,
#        legend.position = "none")
ggsave("../doc/RSPO3_C.pdf",height = 5, width = 5)


ld = read.table("../dat/RSPO3_eBMD/RSPO3.ld")
bim = read.table("../dat/RSPO3_eBMD/RSPO3_UKB.bim")
rownames(bim) = bim$V2
bim = bim[gwas$SNP,]
idxsnp = unlist(lapply(strsplit(cs$cs,'/'),function(x){x[1]}))[4]

pltdat = data.frame(Pos = rep(bim$V4,2),
                     p = c(gwas$mlogP,qtl$mlogP),
                     study = rep(c("BMD GWAS","RSPO3 pQTL"),each=nrow(gwas)),
                     lead = ifelse(gwas$SNP==idxsnp,"Lead SNP","Others"),
                     r2 = rep((ld[,which(gwas$SNP==idxsnp)])^2,2))

ggplot(pltdat, aes(x = Pos / 1000000, y = p, color = r2)) +
  geom_point() +
  facet_grid(study~., scale = "free") +
  scale_color_stepsn(breaks = c(0,0.2,0.4,0.6,0.8,1), colors = c("darkblue","royalblue","green","orange","red")) +
  theme_bw() +
  xlab("Chr6 (Mb)") + 
  #  xlim(127.25,127.75) +
  ylab(expression(-log[10](p-value))) +
  #  scale_shape_manual(values = c(8,16)) +
  labs(size = "", color = expression(r^2)) +
  scale_size_manual(values = c(5,2)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13))
ggsave("../doc/RSPO3_A.pdf",height = 10, width = 7)


ptdat = data.frame(SNP = gwas$SNP, gwas = gwas$mlogP, qtl = qtl$mlogP, r2 = (ld[,which(gwas$SNP==idxsnp)])^2)

ggplot(ptdat, aes(x = qtl, y = gwas, color = r2)) +
  geom_vline(xintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size=2) + theme_bw() +
  scale_color_stepsn(breaks = c(0,0.2,0.4,0.6,0.8,1), colors = c("darkblue","royalblue","green","orange","red")) +
  #scale_color_manual(values = c("hotpink", "darkmagenta", "cyan", "cyan4", "coral", "coral4", "grey")) +
  xlab(expression (-log[10](p-value)~"in"~RSPO3~pQTL)) +
  ylab(expression (-log[10](p-value)~"in"~BMD~GWAS)) +
  labs(color = '') +
  scale_size_manual(values = c(5, 2)) +
  theme(axis.title = element_text (size = 14),
        axis.text = element_text (size = 13),
        legend.text = element_text (size = 13),
        legend.position = "none")
ggsave("../doc/RSPO3_B.pdf",height = 5, width = 5)


## coloc results
coloc = read.table("../dat/RSPO3_eBMD/csusie.coloc.txt", header=T)
coloc_H4 = max(coloc['PP.H4.abf'])
csusie = read.table("../dat/RSPO3_eBMD/csusie.coloc_susie.txt", header=T)
csu_H4 = max(csusie['PP.H4.abf'])
sh_H4 = max(cs$share)
pwcoco = read.table("../dat/RSPO3_eBMD/pwcoco_out.coloc", header=T)
pw_H4 = max(pwcoco$H4)
eca = read.table("../dat/RSPO3_eBMD/BMD_RSPO3_eCAVIAR_col", header=T)
eca_H4 = max(eca$CLPP)

H4dat = data.frame(H4 = c(coloc_H4, csu_H4, eca_H4, pw_H4, sh_H4), 
                   method = c("COLOC","COLOC+SuSiE","eCAVIAR","PWCoCo","SharePro"))

ggplot(H4dat,aes(x=method,y=H4,shape=method)) + 
  geom_point(size = 4) + 
  theme_bw() + 
  xlab("") +
  ylab("Colocalization probability") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  scale_y_continuous(breaks = c(0,0.2,0.8,1)) + 
  scale_shape_manual(values = c(15,17,18,19,3)) +
  labs(shape="")
ggsave("../doc/RSPO3_D.pdf",height = 3, width = 10)

