## Run GWAS 
library(dplyr)
library(bigsnpr)
library(readr)
library(mashr)
library(AnnotationDbi)
library(optparse)
library(collections)
library(bigparallelr)
library('matrixcalc')
assert_cores(1)

option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="which site to use",metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}

pheno = opt$pheno

betas = read.table(paste0('combined_mash_effects/',pheno,'.all.betas.strong.txt'))
stderrs = read.table(paste0('combined_mash_stderrs/',pheno,'.all.stderrs.strong.txt'))

mat_betas = matrix(unlist(betas),ncol=8,byrow = TRUE)
mat_stderrs = matrix(unlist(stderrs),ncol=8,byrow = TRUE)

sites = c('MI','MO','NE','OK','SD','TX1','TX2','TX3')
data=mash_set_data(mat_betas,mat_stderrs) 

canonical = cov_canonical(data = data)
matnames = ls(canonical)
for (matrix in matnames) {
  if (matrix %in% c('equal_effects','identity','simple_het_1','simple_het_2','simple_het_3')) {
  matrixout = matrix(unlist(canonical[matrix]),ncol=8,byrow=TRUE)
  row.names(matrixout) = sites
  colnames(matrixout) = sites
  write.csv(matrixout, file = paste0('matrices/',pheno,'.',matrix,'.csv'))
  }
}
print(canonical)
for (i in 1:8) {
  singleton = paste0('singletons_',i)
  print(singleton)
  singmat = canonical[singleton]
  site = sites[i]
  
  singmat = matrix(unlist(canonical[singleton]),ncol=8,byrow=TRUE)
  row.names(singmat) = sites
  colnames(singmat) = sites
  write.csv(singmat, file = paste0('matrices/',pheno,'.singleton.',site,'.csv'))
}
# write.table(c(likelihood),file = paste0('likelihoods/step',step,'/',pheno,'.',outname,'.',cohort,'.',effect,'.txt'),row.names = FALSE,col.names = FALSE)
# write.table(weights,file = paste0('likelihoods/step',step,'/',pheno,'.',outname,'.',cohort,'.',effect,'.weights.txt'),col.names=FALSE)
