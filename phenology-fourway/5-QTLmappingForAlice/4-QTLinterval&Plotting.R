setwd("~/Github/pvdiv-phenology-gxe/phenology-fourway/5-QTLmappingForAlice/")
library(qtl)
library(qtl2)
library(plyr)
library(tidyverse)
library(qtl2convert)

############# WORK ON GENO
rm(list=ls())
load('cross_reduced_marker_FL50_phenotypes_Alice_7_SITES.RData')
load('scan1Results_FL50_phenotypes_Alice_7_SITES.RData')

map  =  lapply(pull.map(cross), function(x) x[1,])
###get the threshold value for each model and output the QTL intervals
threshold_full = summary(s1permfull,alpha = 0.05)
threshold_reduced = summary(s1permreduced,alpha = 0.05)

###  g x e
gxe = s1full - s1reduced
gxe_perm = threshold_full - threshold_reduced ###or gxe_perm = summary(gxe, alpha=0.05)

###get the peak for full model and gxe model
peaks_full = find_peaks(s1full, map, threshold = threshold_full, drop = 1.5, peakdrop = 5)
peaks_gxe = find_peaks(gxe, map, threshold = gxe_perm, drop = 1.5)

flank_marker_full = peaks_full %>% mutate(flank_lo = find_marker(map, chr, ci_lo), flank_hi=find_marker(map, chr, ci_hi))
write.csv(flank_marker_full,'QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_Full_model.csv') ###change the name you want

flank_marker_gxe = peaks_gxe %>% mutate(flank_lo = find_marker(map, chr, ci_lo), flank_hi=find_marker(map, chr, ci_hi))
write.csv(flank_marker_gxe,'QTLsWithFlankMarkers_FL50_phenotypes_Alice_7_SITES_gxe_model.csv')

###########
png('../2-figures/FL50_7_SITES_LOD_FL50_phenotypes_Alice_7_SITES.png', width = 800)
plot(s1full, map, lodcolumn=2, col="slateblue", ylim=c(0,56))
plot(s1full,map, lodcolumn=3, col="violetred", add=T)
plot(s1full,map, lodcolumn=4, col="darkgreen", add=T)
plot(s1full,map, lodcolumn=5, col="yellow", add=T)
abline(h=threshold_full[2],lty=2,lwd=3,col=c("slateblue"))
abline(h=threshold_full[3],lty=2,lwd=3,col=c("violetred"))
abline(h=threshold_full[4],lty=2,lwd=3,col=c("darkgreen"))
abline(h=threshold_full[5],lty=2,lwd=3,col=c("yellow"))

legend("topleft", lwd=2, col=c("slateblue", "violetred","darkgreen","yellow"), colnames(s1full)[c(2,3,4,5)], bg="gray90")
dev.off()


flank_marker_full = flank_marker_full %>% mutate(GXE= ifelse(pos %in% flank_marker_gxe$pos, 'Y','N'))%>%
  mutate(ChrNo = str_extract(chr, '[1-9]'))%>% mutate_at(vars(ChrNo),as.numeric) %>%
  mutate(Dist = case_when(
    str_detect(chr,'K') ~ (ChrNo*2-1),
    str_detect(chr,'N') ~ ChrNo*2
  ))
# ###6. plot QTL on genetic map
# load('sexAveraged.cross.rda')
png('../2-figures/QTLs_FL50_phenotypes_Alice_7_SITES.png',width = 900)
cols  =  c("violetred","yellow","darkgreen","slateblue") ###add your colors depending on how many traits you have
library(qtlTools)
with(flank_marker_full, segmentsOnMap(cross, phe = lodcolumn, chr = chr,
                                      l = ci_lo, h = ci_hi,
                                      peaklod = lod, peakcM = pos,
                                      showPeaks = TRUE,
                                      chrBuffer = c(0.15,0.1),
                                      tick.width=0.05, lwd="byLod",
                                      col = cols,
                                      leg.inset=.57, legendCex=1,
                                      legendPosition="topleft")
)
#mtext('QTL_7_SITES',side = 3)

for (i in 1:nrow(flank_marker_full)){
    if (flank_marker_full$GXE[i]=='Y')
    text(x=(flank_marker_full$Dist[i]+0.4), y=(flank_marker_full$ci_lo[i]),'*',col = 'red',cex = 1.5)
  }

dev.off()

