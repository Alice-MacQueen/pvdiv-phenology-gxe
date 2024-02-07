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
  make_option(c("-p","--pheno"), type="character", default=NULL,help="which site to use",metavar="character"),
  make_option(c("-c","--cohort"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-e","--effects"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}

pheno = opt$pheno
cohort = opt$cohort
effect = opt$effects

betas = read.table(paste0('combined_mash_effects/',pheno,'.',cohort,'.betas.',effect,'.txt'))
stderrs = read.table(paste0('combined_mash_stderrs/',pheno,'.',cohort,'.stderrs.',effect,'.txt'))

mat_betas = matrix(unlist(betas),ncol=8,byrow = TRUE)
mat_stderrs = matrix(unlist(stderrs),ncol=8,byrow = TRUE)

data=mash_set_data(mat_betas,mat_stderrs) 

canonical = cov_pca(data,5)
sites = c('MI','MO','NE','OK','SD','TX1','TX2','TX3')

matnames = ls(canonical)
for (mat in matnames) {
  matrixout = matrix(unlist(canonical[mat]),ncol=8,byrow=TRUE)
  row.names(matrixout) = sites
  colnames(matrixout) = sites
  print(matrixout)
  write.csv(matrixout, file = paste0('matrices/',pheno,'.',cohort,'.',effect,'.',mat,'.csv'))
}
