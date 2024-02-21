library(ggplot2)
setwd('~/Desktop/Projects/SharePro_loc/plt')

res = data.frame(trait=c('LDL', 'LDL', 'Tg', 'Tg', 'Tg'), 
                 protein=c('PCSK9', 'APOB', 'ANGPTL3', 'LPL', 'APOC3'))

res$pair = paste0(res$trait,':',res$protein)

coloc_H4 = c()
csu_H4 = c()
eca_H4 = c()
pw_H4 = c()
sh_H4 = c()
prs = c('1e-3', '1e-4', '1e-5', '1e-6', '1e-7')
for (pr in prs) {
  for (i in 1:nrow(res)) {
    pwcoco = read.table(paste0('../dat/LDL/res_', pr, '/pwcoco_', res$trait[i], '_', res$protein[i], '.coloc'), header=T)
    coloc_H4 = c(coloc_H4, pwcoco[1, 'H4'])
    pw_H4 = c(pw_H4, max(pwcoco$H4))
    csu_file = paste0('../dat/LDL/res_', pr, '/csusie_', res$trait[i], '_', res$protein[i], '.txt')
    if (file.exists(csu_file)) {
      csusie = read.table(csu_file, header=T)
      csu_H4 = c(csu_H4, max(csusie$PP.H4.abf))
    } else {
      csu_H4 = c(csu_H4, 0)
    }
    eca = read.table(paste0('../dat/LDL/res_1e-5/ecaviar_', res$trait[i], '_', res$protein[i], '_col'), header=T)
    eca_H4 = c(eca_H4, max(eca$CLPP))
    sh = read.table(paste0('../dat/LDL/res_', pr, '/SH_', res$trait[i], '.', res$protein[i], '.txt'), header=T)
    sh_H4 = c(sh_H4, max(sh$share))
  }
}

# Supplementary Table S5
df_res = data.frame(prior = rep(prs, each=5), protein=rep(res$pair, 5), coloc = coloc_H4, csusie = csu_H4, ecaviar = eca_H4, pwcoco = pw_H4, sharepro = sh_H4)
write.table(df_res,"../doc/sharepro_loc_ldl.csv",sep=',', row.names = F, quote = F)


pltdat = data.frame(prior = rep(rep(prs, each=5), 5), 
                    protein=rep(rep(res$pair, 5), 5),
                    method = rep(c('COLOC', 'COLOC+SuSiE', 'eCAVIAR', 'PWCoCo', 'SharePro'), each=5*5),
                    H4 = c(coloc_H4, csu_H4, eca_H4, pw_H4, sh_H4))
pltdat$prior = factor(pltdat$prior, levels = c('1e-7','1e-6','1e-5','1e-4','1e-3'))
pltdat$protein = factor(pltdat$protein, levels=c('LDL:APOB', 'LDL:PCSK9', 'Tg:LPL', 'Tg:APOC3', 'Tg:ANGPTL3'))
ggplot(pltdat, aes(x=prior, y=H4, shape=method)) + 
  geom_point(position = position_dodge(width = 0.5), size = 3) + theme_bw() +
  facet_grid(protein~.) + 
  scale_y_continuous(breaks = c(0,0.2,0.8,1)) +
  scale_x_discrete(labels = c(expression(10^-7),expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3))) +
  scale_shape_manual(values = c(15,17,18,19,3)) +
  ylab('Colocalization probability') +
  labs(shape = "") + 
  theme(axis.title = element_text(size = 19),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        title = element_text(size = 14))
ggsave("../doc/LDL.pdf",height = 10, width = 10)


