library(ggplot2)
setwd('~/Desktop/Projects/SharePro_loc/plt')
colors = c("hotpink", "darkmagenta", "cyan", "cyan4", "coral", "coral4", "plum1", "purple", "olivedrab", "seagreen")

#args = commandArgs(trailingOnly=T)
trait = 'Tg'
protein = 'LPL'
qdir = paste0('../dat/LDL/ss/', trait, '.', protein, '.', protein, '.pwcoco')
gdir = paste0('../dat/LDL/ss/', trait, '.', protein, '.', trait, '.pwcoco')
lddir = paste0('../dat/LDL/ss/', protein, '.ld')
csdir = paste0('../dat/LDL/res_1e-5/SH_', trait, '.', protein, '.txt')
bimdir = paste0('../dat/LDL/bed/', protein, '.bim')

gwas = read.table(gdir, header=T)
qtl = read.table(qdir, header=T)
cs = read.table(csdir, header=T)
gwas$mlogP = -log10(exp(1)) * pchisq((gwas$BETA/gwas$SE)^2, df=1, lower.tail = F, log.p=T)
qtl$mlogP = -log10(exp(1)) * pchisq((qtl$BETA/qtl$SE)^2, df=1, lower.tail = F, log.p=T)

plotdat = data.frame(SNP = gwas$SNP, gwas = gwas$mlogP, qtl = qtl$mlogP, cs = "others")
for (i in 1:nrow(cs)) {
  cslist <- strsplit(cs$cs[i], '/')[[1]]
  plotdat$cs[plotdat$SNP %in% cslist] <- paste("Effect group", i)
}
plotdat$cs = factor(plotdat$cs, levels = c(paste("Effect group", c(1:nrow(cs))), "others"))
plotdat = plotdat[order(plotdat$cs, decreasing = T),]
ggplot(plotdat, aes(x = qtl, y = gwas, color = cs)) +
  geom_vline(xintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point() + theme_bw() +
  scale_color_manual(values = c(colors[1:nrow(cs)], 'grey')) +
  xlab(expression(-log[10](p-value)~"in"~pQTL)) +
  ylab(expression(-log[10](p-value)~"in"~GWAS)) +
  labs(color = '') +
  scale_size_manual(values = c(5, 2)) +
  theme(axis.title = element_text (size = 14),
        axis.text = element_text (size = 13),
        legend.text = element_text (size = 13),
        legend.position = "none") -> plt
ggsave(paste0("../doc/", protein, "_", trait, "_C.pdf"), plt, height = 5, width = 5)


ld = read.table(lddir)
bim = read.table(bimdir)
rownames(bim) = bim$V2
bim = bim[gwas$SNP,]
idxsnp = unlist(lapply(strsplit(cs$cs,'/'),function(x){x[1]}))[1]


pltdat = data.frame(Pos = rep(bim$V4,2),
                    p = c(gwas$mlogP,qtl$mlogP),
                    study = rep(c(paste0(trait, " GWAS"), paste0(protein, " pQTL")),each=nrow(gwas)),
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
        strip.text = element_text(size = 13),
        legend.position = "none") -> plt
ggsave(paste0("../doc/", protein, "_", trait, "_A.pdf"), plt, height = 7, width = 4)



ptdat = data.frame(SNP = gwas$SNP, gwas = gwas$mlogP, qtl = qtl$mlogP, r2 = (ld[,which(gwas$SNP==idxsnp)])^2)

ggplot(ptdat, aes(x = qtl, y = gwas, color = r2)) +
  geom_vline(xintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size=2) + theme_bw() +
  scale_color_stepsn(breaks = c(0,0.2,0.4,0.6,0.8,1), colors = c("darkblue","royalblue","green","orange","red")) +
  #scale_color_manual(values = c("hotpink", "darkmagenta", "cyan", "cyan4", "coral", "coral4", "grey")) +
  xlab(expression (-log[10](p-value)~"in"~pQTL)) +
  ylab(expression (-log[10](p-value)~"in"~GWAS)) +
  labs(color = '') +
  scale_size_manual(values = c(5, 2)) +
  theme(axis.title = element_text (size = 14),
        axis.text = element_text (size = 13),
        legend.text = element_text (size = 13),
        legend.position = "none") -> plt
ggsave(paste0("../doc/", protein, "_", trait, "_B.pdf"), plt, height = 7, width = 7)
