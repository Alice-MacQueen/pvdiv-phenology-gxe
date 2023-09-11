## Run GWAS 
library(dplyr)
library(bigsnpr)
library(readr)
library(mashr)
library(AnnotationDbi)
library(optparse)
library(collections)
library(bigparallelr)

assert_cores(1)

option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="which site to use",metavar="character"),
  make_option(c("-c","--cohort"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-m","--matrices"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-s","--step"), type="character", default=NULL,help="step in the greedy algorithm",metavar="character"),
  make_option(c("-e","--effects"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-b","--bootstrap"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}

pheno = opt$pheno
cohort = opt$cohort
matrices = unlist(strsplit(opt$matrices,','))
effect = opt$effects
step = opt$step
bootstrap=opt$bootstrap
print(matrices)
betas = read.table(paste0('bootstrap_combined_mash_effects/',pheno,'.',cohort,'.betas.',effect,'.',bootstrap,'.txt'))
stderrs = read.table(paste0('bootstrap_combined_mash_stderrs/',pheno,'.',cohort,'.stderrs.',effect,'.',bootstrap,'.txt'))

mat_betas = matrix(unlist(betas),ncol=8,byrow = TRUE)
mat_stderrs = matrix(unlist(stderrs),ncol=8,byrow = TRUE)

data=mash_set_data(mat_betas,mat_stderrs) 


covmats = list()
labels = list()
for (matrix in matrices) {
  label = unlist(strsplit(toString(matrix),'.',fixed = TRUE))
  label=paste(label[2],label[3],sep='.')
  temp = read.csv(matrix,header = TRUE, row.names = 1)
  temp = matrix(unlist(temp),ncol=8,byrow = TRUE)
  covmats[[label]] = temp
  labels = append(labels,label)
}

print(typeof(matrices))
g = mash(data,covmats)
likelihood = mash_compute_loglik(g,data)
print(get_estimated_pi(g))
outname = tail(labels,n=1)
write.table(c(likelihood),file = paste0('bootstrap_likelihoods/step',step,'/',pheno,'.',outname,'.',cohort,'.',effect,'.',bootstrap,'.txt'),row.names = FALSE,col.names = FALSE)

