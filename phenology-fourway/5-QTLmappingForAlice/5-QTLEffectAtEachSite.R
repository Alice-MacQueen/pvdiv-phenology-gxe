###QTL Effect plot
rm(list=ls())
library(qtl)
library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("~/Github/pvdiv-phenology-gxe/phenology-fourway/5-QTLmappingForAlice/")
###load the original cross and reduce the markers from 4700 to 1185
load('sexAveraged.cross.rda', verbose = T)
map  =  lapply(pull.map(cross), function(x) x[1,])
# reduced_markers  =  unlist(
#   lapply(
#     reduce_markers(
#       lapply(pull.map(cross), function(x)
#         x[1,]), min_distance = 1), names))
#
# mars2drop  =  markernames(cross)[!markernames(cross) %in% reduced_markers]
# cross  =  drop.markers(cross, mars2drop)
cross = calc.genoprob(cross, step=1, error.prob = 1e-3,map.function = "kosambi",stepwidth = "max")
######### Converting probalibity from qtl1 to qtl2
probs  =  probs_qtl_to_qtl2(cross)$probs
##calculate the kinship matrix for single site
kinship_loco_Single = calc_kinship(probs,use_allele_probs = F,type = "loco")

###extract the genotype ID from the original cross file, so we can align the phenology file with genotype ID
ID = data.frame(ID=cross$pheno$id)
ID = ID %>% mutate_at(vars('ID'), list(as.character)) %>% mutate_at(vars('ID'), list(as.numeric))

######get the phenotype at each site
phenology = read.csv('FL50_phenotypes_Alice.csv')
colnames(phenology)[1:2] = c('ID','Covariate')
str(phenology) ###make sure your phenology traits are numeric, if not, change them into numeric
phenos = phenology %>% mutate_at(vars('ID'), list(as.character)) %>%
  mutate_at(vars('ID'), list(as.numeric))%>% group_by(Covariate,ID ) %>%
  summarize_all(mean, na.rm=T) %>% ungroup()%>% as.data.frame()

###QTL effect for each trait---FL50
EffectMatrix = vector('list', length=0)
t=0
for (envi in unique(phenos$Covariate)){
    t=t+1
    pheno = phenos %>% filter(Covariate==envi) %>% dplyr::select(-Covariate) %>%
      select(ID, FL50) %>% right_join(ID) %>%  column_to_rownames('ID') %>% as.matrix()

    ##ran effect scan for each chromosome
    eff_out = vector('list',length = 0)
    se_out = vector('list',length = 0)
    add_contrasts=cbind(mu=c(1,1,1,1), add1=c(-1, 1, -1, 1),add2=c(-1, -1, 1, 1), dom=c(-1, 1, 1, -1))

    z = 0
    for (chr in chrnames(cross)){
      print(c(chr, envi,t))
      z= z+1
      eff = scan1coef(probs[,chr],pheno,kinship = kinship_loco_Single[chr], se=T, contrasts = add_contrasts)
      se = attr(eff, 'SE')
      eff_out[[z]] =  eff %>% as.data.frame() %>% rownames_to_column('MARKER')
      se_out[[z]] =  se %>% as.data.frame() %>% rownames_to_column('MARKER')
    }
    eff_out = bind_rows(eff_out)
    se_out = bind_rows(se_out)

    se_out2 = se_out %>% dplyr::select(-mu) %>% gather('CROSS','SE',-MARKER)
    eff_out2 = eff_out %>% dplyr::select(-mu) %>% gather('CROSS','EFF',-MARKER) %>% mutate(SE=se_out2$SE, Covariate=envi)

    EffectMatrix[[t]] = eff_out2
  }

EffectMatrix2 = bind_rows(EffectMatrix)

###unlist nmap and merge with Effect file
nmap = data.frame(pos=unlist(map))
nnmap = nmap %>% rownames_to_column('MARKER') %>% rowwise() %>%
  mutate(chr = str_sub(MARKER, 1,2), MARKER=str_sub(MARKER, 4,nchar(MARKER)))%>%
  ungroup() %>% mutate(pos=round(pos,4))

flank_marker_full = read.csv('QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_Full_model.csv')
flank_marker_full = flank_marker_full %>% mutate_at(vars(chr),list(as.character))

for (i in 1:nrow(flank_marker_full)){
  tmp = nnmap %>% filter(chr==flank_marker_full$chr[i] & pos==round(flank_marker_full$pos[i],4))
  flank_marker_full$MARKER[i] = tmp$MARKER
}

QTL_full_model = flank_marker_full %>% dplyr::rename(trait=lodcolumn)

QTL_Eff = QTL_full_model %>% left_join(EffectMatrix2)%>%
  mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>%
  mutate(CROSS = str_replace(CROSS,'add1','A x B'))


QTL_Eff$Covariate = factor(QTL_Eff$Covariate, levels=c('KING','PKLE','STIL','CLMB',
                                                       'LINC','KBSM','BRKG'))


ggplot(data=(QTL_Eff %>% filter(trait=='FL50')), aes(Covariate, EFF, group=CROSS)) + geom_point(size=2) + facet_grid(CROSS~MARKER) +
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.6,width=0.5) +
  ylab('QTL Effect') +xlab("") +coord_flip()+theme_bw()+
  theme(text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))

ggsave('../2-figures/QTL_Effects_FL50_with_dom_2019.png', width = 11)

###QTL effect for each trait---CGDD_12C
EffectMatrix = vector('list', length=0)
t=0
for (envi in unique(phenos$Covariate)){
  t=t+1
  pheno = phenos %>% filter(Covariate==envi) %>% dplyr::select(-Covariate) %>%
    select(ID, CGDD_12C) %>% right_join(ID) %>%  column_to_rownames('ID') %>% as.matrix()

  ##ran effect scan for each chromosome
  eff_out = vector('list',length = 0)
  se_out = vector('list',length = 0)
  add_contrasts=cbind(mu=c(1,1,1,1), add1=c(-1, 1, -1, 1),add2=c(-1, -1, 1, 1), dom=c(-1, 1, 1, -1))

  z = 0
  for (chr in chrnames(cross)){
    print(c(chr, envi,t))
    z= z+1
    eff = scan1coef(probs[,chr],pheno,kinship = kinship_loco_Single[chr], se=T, contrasts = add_contrasts)
    se = attr(eff, 'SE')
    eff_out[[z]] =  eff %>% as.data.frame() %>% rownames_to_column('MARKER')
    se_out[[z]] =  se %>% as.data.frame() %>% rownames_to_column('MARKER')
  }
  eff_out = bind_rows(eff_out)
  se_out = bind_rows(se_out)

  se_out2 = se_out %>% dplyr::select(-mu) %>% gather('CROSS','SE',-MARKER)
  eff_out2 = eff_out %>% dplyr::select(-mu) %>% gather('CROSS','EFF',-MARKER) %>% mutate(SE=se_out2$SE, Covariate=envi)

  EffectMatrix[[t]] = eff_out2
}

EffectMatrix2 = bind_rows(EffectMatrix)

###unlist nmap and merge with Effect file
nmap = data.frame(pos=unlist(map))
nnmap = nmap %>% rownames_to_column('MARKER') %>% rowwise() %>%
  mutate(chr = str_sub(MARKER, 1,2), MARKER=str_sub(MARKER, 4,nchar(MARKER)))%>%
  ungroup() %>% mutate(pos=round(pos,4))

flank_marker_full = read.csv('QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_Full_model.csv')
flank_marker_full = flank_marker_full %>% mutate_at(vars(chr),list(as.character))

for (i in 1:nrow(flank_marker_full)){
  tmp = nnmap %>% filter(chr==flank_marker_full$chr[i] & pos==round(flank_marker_full$pos[i],4))
  flank_marker_full$MARKER[i] = tmp$MARKER
}

QTL_full_model = flank_marker_full %>% dplyr::rename(trait=lodcolumn)

QTL_Eff = QTL_full_model %>% left_join(EffectMatrix2)%>%
  mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>%
  mutate(CROSS = str_replace(CROSS,'add1','A x B'))


QTL_Eff$Covariate = factor(QTL_Eff$Covariate, levels=c('KING','PKLE','STIL','CLMB',
                                                       'LINC','KBSM','BRKG'))


ggplot(data=(QTL_Eff %>% filter(trait=='CGDD_12C')), aes(Covariate, EFF, group=CROSS)) + geom_point(size=2) + facet_grid(CROSS~MARKER) +
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.6,width=0.5) +
  ylab('QTL Effect') +xlab("") +coord_flip()+theme_bw()+
  theme(text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))

ggsave('../2-figures/QTL_Effects_CGDD_12C_with_Dom_2019.png', width = 11)



###QTL effect for each trait---dyln_fl50
EffectMatrix = vector('list', length=0)
t=0
for (envi in unique(phenos$Covariate)){
  t=t+1
  pheno = phenos %>% filter(Covariate==envi) %>% dplyr::select(-Covariate) %>%
    select(ID, dyln_fl50) %>% right_join(ID) %>%  column_to_rownames('ID') %>% as.matrix()

  ##ran effect scan for each chromosome
  eff_out = vector('list',length = 0)
  se_out = vector('list',length = 0)
  add_contrasts=cbind(mu=c(1,1,1,1), add1=c(-1, 1, -1, 1),add2=c(-1, -1, 1, 1), dom=c(-1, 1, 1, -1))

  z = 0
  for (chr in chrnames(cross)){
    print(c(chr, envi,t))
    z= z+1
    eff = scan1coef(probs[,chr],pheno,kinship = kinship_loco_Single[chr], se=T, contrasts = add_contrasts)
    se = attr(eff, 'SE')
    eff_out[[z]] =  eff %>% as.data.frame() %>% rownames_to_column('MARKER')
    se_out[[z]] =  se %>% as.data.frame() %>% rownames_to_column('MARKER')
  }
  eff_out = bind_rows(eff_out)
  se_out = bind_rows(se_out)

  se_out2 = se_out %>% dplyr::select(-mu) %>% gather('CROSS','SE',-MARKER)
  eff_out2 = eff_out %>% dplyr::select(-mu) %>% gather('CROSS','EFF',-MARKER) %>% mutate(SE=se_out2$SE, Covariate=envi)

  EffectMatrix[[t]] = eff_out2
}

EffectMatrix2 = bind_rows(EffectMatrix)

###unlist nmap and merge with Effect file
nmap = data.frame(pos=unlist(map))
nnmap = nmap %>% rownames_to_column('MARKER') %>% rowwise() %>%
  mutate(chr = str_sub(MARKER, 1,2), MARKER=str_sub(MARKER, 4,nchar(MARKER)))%>%
  ungroup() %>% mutate(pos=round(pos,4))

flank_marker_full = read.csv('QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_Full_model.csv')
flank_marker_full = flank_marker_full %>% mutate_at(vars(chr),list(as.character))

for (i in 1:nrow(flank_marker_full)){
  tmp = nnmap %>% filter(chr==flank_marker_full$chr[i] & pos==round(flank_marker_full$pos[i],4))
  flank_marker_full$MARKER[i] = tmp$MARKER
}

QTL_full_model = flank_marker_full %>% dplyr::rename(trait=lodcolumn)

QTL_Eff = QTL_full_model %>% left_join(EffectMatrix2)%>%
  mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>%
  mutate(CROSS = str_replace(CROSS,'add1','A x B'))


QTL_Eff$Covariate = factor(QTL_Eff$Covariate, levels=c('KING','PKLE','STIL','CLMB',
                                                       'LINC','KBSM','BRKG'))


ggplot(data=(QTL_Eff %>% filter(trait=='dyln_fl50')), aes(Covariate, EFF, group=CROSS)) + geom_point(size=2) + facet_grid(CROSS~MARKER) +
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.6,width=0.5) +
  ylab('QTL Effect') +xlab("") +coord_flip()+theme_bw()+
  theme(text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))

ggsave('../2-figures/QTL_Effects_dyln_fl50_with_Dom_2019.png', width = 11)


###QTL effect for each trait---dyln_change_sec
EffectMatrix = vector('list', length=0)
t=0
for (envi in unique(phenos$Covariate)){
  t=t+1
  pheno = phenos %>% filter(Covariate==envi) %>% dplyr::select(-Covariate) %>%
    select(ID, dyln_change_sec) %>% right_join(ID) %>%  column_to_rownames('ID') %>% as.matrix()

  ##ran effect scan for each chromosome
  eff_out = vector('list',length = 0)
  se_out = vector('list',length = 0)
  add_contrasts=cbind(mu=c(1,1,1,1), add1=c(-1, 1, -1, 1),add2=c(-1, -1, 1, 1), dom=c(-1, 1, 1, -1))

  z = 0
  for (chr in chrnames(cross)){
    print(c(chr, envi,t))
    z= z+1
    eff = scan1coef(probs[,chr],pheno,kinship = kinship_loco_Single[chr], se=T, contrasts = add_contrasts)
    se = attr(eff, 'SE')
    eff_out[[z]] =  eff %>% as.data.frame() %>% rownames_to_column('MARKER')
    se_out[[z]] =  se %>% as.data.frame() %>% rownames_to_column('MARKER')
  }
  eff_out = bind_rows(eff_out)
  se_out = bind_rows(se_out)

  se_out2 = se_out %>% dplyr::select(-mu) %>% gather('CROSS','SE',-MARKER)
  eff_out2 = eff_out %>% dplyr::select(-mu) %>% gather('CROSS','EFF',-MARKER) %>% mutate(SE=se_out2$SE, Covariate=envi)

  EffectMatrix[[t]] = eff_out2
}

EffectMatrix2 = bind_rows(EffectMatrix)

###unlist nmap and merge with Effect file
nmap = data.frame(pos=unlist(map))
nnmap = nmap %>% rownames_to_column('MARKER') %>% rowwise() %>%
  mutate(chr = str_sub(MARKER, 1,2), MARKER=str_sub(MARKER, 4,nchar(MARKER)))%>%
  ungroup() %>% mutate(pos=round(pos,4))

flank_marker_full = read.csv('QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_Full_model.csv')
flank_marker_full = flank_marker_full %>% mutate_at(vars(chr),list(as.character))

for (i in 1:nrow(flank_marker_full)){
  tmp = nnmap %>% filter(chr==flank_marker_full$chr[i] & pos==round(flank_marker_full$pos[i],4))
  flank_marker_full$MARKER[i] = tmp$MARKER
}

QTL_full_model = flank_marker_full %>% dplyr::rename(trait=lodcolumn)

QTL_Eff = QTL_full_model %>% left_join(EffectMatrix2)%>%
  mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>%
  mutate(CROSS = str_replace(CROSS,'add1','A x B'))


QTL_Eff$Covariate = factor(QTL_Eff$Covariate, levels=c('KING','PKLE','STIL','CLMB',
                                                       'LINC','KBSM','BRKG'))


ggplot(data=(QTL_Eff %>% filter(trait=='dyln_change_sec')), aes(Covariate, EFF, group=CROSS)) + geom_point(size=2) + facet_grid(CROSS~MARKER) +
  geom_hline(yintercept = 0,linetype=2)+ xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.6,width=0.5) +
  ylab('QTL Effect') +xlab("") +coord_flip()+theme_bw()+
  theme(text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))

ggsave('../2-figures/QTL_Effects_dyln_change_sec_with_Dom_2019.png', width = 11)

#####End




