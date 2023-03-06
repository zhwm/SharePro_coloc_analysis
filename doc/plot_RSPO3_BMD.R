setwd("~/Desktop/GitHub/SharePro_coloc_analysis/dat/")
library(ggplot2)

gwas <- read.table("BMD_SH.txt")
gwas$mlogP <- -log10(exp(1)) * pchisq(gwas$V2^2,df=1,lower.tail = F,log.p = T)
qtl <- read.table("RSPO3_SH.txt")
qtl$mlogP <- -log10(exp(1)) * pchisq(qtl$V2^2,df=1,lower.tail = F,log.p = T)
ld <- read.table("RSPO3.ld")
bim <- read.table("RSPO3_UKB.bim")
rownames(bim) <- bim$V2
bim <- bim[gwas$V1,]
cs <- read.table("BMD_SH.txt_RSPO3_SH.txt.cs",header=T)
cs$topsnp <- unlist(lapply(strsplit(cs$cs,'/'),function(x){x[1]}))
colnames(gwas)[1] <- "SNP"

pltdat <- data.frame(Pos = rep(bim$V4,2),
                     p = c(gwas$mlogP,qtl$mlogP),
                     study = rep(c("eBMD GWAS","RSPO3 pQTL"),each=nrow(gwas)),
                     lead = ifelse(gwas$SNP==cs$topsnp[2],"Lead SNP","Others"),
                     r2 = rep((ld[,which(gwas$SNP==cs$topsnp[2])])^2,2))

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
        strip.text = element_text(size = 13)) -> plt
ggsave("~/Desktop/GitHub/SharePro_coloc_analysis/doc/BMD_RSPO3.pdf",plt,height = 10, width = 10)

H4dat <- data.frame(H4=c(0.000224131,
                         2.95e-17,
                         0.2040571,
                         9.53551e-12,
                         1.0),
                    method=c("eCAVIAR",
                             "COLOC",
                             "COLOC+SuSiE",
                             "PWCoCo",
                             "SharePro"))

ggplot(H4dat,aes(x=method,y=H4)) + 
  geom_point(size = 4) + 
  theme_bw() + 
  xlab("") +
  ylab("Colocalization probability") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  scale_y_continuous(breaks = c(0,0.2,0.8,1))
ggsave("~/Desktop/GitHub/SharePro_coloc_analysis/doc/BMD_RSPO3_H4.pdf", height = 3, width = 10)

## Format supplementary tables
# SharePro
prior <- c('1e-3/','1e-4/','1e-5/','1e-6/','1e-7/')
ddir <- '~/Desktop/GitHub/SharePro_coloc_analysis/dat/sensitivity/'

for (i in prior) {
  cs <- read.table(paste0(ddir,i,"BMD_SH.txt_RSPO3_SH.txt.cs"),header=T)
  cs$topsnp <- unlist(lapply(strsplit(cs$cs,'/'),function(x){x[1]}))
  cs$`Effect top variant` = cs$topsnp
  cs$`Colocalization Prob` = cs$share
  cs$`BMD Causal Prob` = unlist(lapply(strsplit(cs$causalProb,','),function(x){x[1]}))
  cs$`RSPO3 Causal Prob` = unlist(lapply(strsplit(cs$causalProb,','),function(x){x[2]}))
  cs$`BMD Specific Prob` = unlist(lapply(strsplit(cs$specific,','),function(x){x[1]}))
  cs$`RSPO3 Specific Prob` = unlist(lapply(strsplit(cs$specific,','),function(x){x[2]}))
  cs$`BMD effect size` = unlist(lapply(strsplit(cs$beta,','),function(x){x[1]}))
  cs$`RSPO3 effect size` = unlist(lapply(strsplit(cs$beta,','),function(x){x[2]}))
  cs$`Effect composition` = cs$cs
  cs$`Effect composition Prob` = cs$variantProb
  write.table(cs[,c("Effect top variant",
                    "BMD Causal Prob","BMD Specific Prob","BMD effect size",
                    "RSPO3 Causal Prob","RSPO3 Specific Prob","RSPO3 effect size",
                    "Colocalization Prob",
                    "Effect composition",
                    "Effect composition Prob")],paste0("~/Desktop/GitHub/SharePro_coloc_analysis/doc/BMD_RSPO3",gsub("/","",i),'.csv'),
              sep=',',col.names = T, row.names = F)
}

library(corrplot)
corrplot(as.matrix(ld^2),
         tl.pos = "n",col.lim = c(0,1),
         col = c("grey","grey","grey","grey","darkblue","royalblue","green","orange","red"),
         addgrid.col = NA)
