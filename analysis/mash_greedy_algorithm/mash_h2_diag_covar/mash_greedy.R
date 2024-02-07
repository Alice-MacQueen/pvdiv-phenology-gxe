## Run GWAS 
suppressMessages(library(dplyr))
suppressMessages(library(bigsnpr))
suppressMessages(library(readr))
suppressMessages(library(mashr))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(optparse))
suppressMessages(library(collections))
suppressMessages(library(bigparallelr))
suppressMessages(library('matrixcalc'))
assert_cores(1)

option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="which site to use",metavar="character"),
  make_option(c("-c","--cohort"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-m","--matrices"), type="character", default=NULL,help="gulf, midwest, or all",metavar="character"),
  make_option(c("-s","--step"), type="character", default=NULL,help="step in the greedy algorithm",metavar="character"),
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
matrices = unlist(strsplit(opt$matrices,','))
effect = opt$effects
step = opt$step

betas = read.table(paste0('combined_mash_effects/',pheno,'.',cohort,'.betas.',effect,'.txt'))
stderrs = read.table(paste0('combined_mash_stderrs/',pheno,'.',cohort,'.stderrs.',effect,'.txt'))

mat_betas = matrix(unlist(betas), ncol = 8, byrow = TRUE)
mat_stderrs = matrix(unlist(stderrs), ncol = 8, byrow = TRUE)

data = mash_set_data(mat_betas, mat_stderrs) 


covmats = list()
labels = list()
for (matrix in matrices) {

  label=gsub('.csv','',matrix)
  label = gsub('matrices/','',label)
  temp = read.csv(matrix,header = TRUE, row.names = 1)
  temp = matrix(unlist(temp),ncol=8,byrow = TRUE)
  covmats[[label]] = temp
  labels = append(labels,label)
  
}

g = mash(data,covmats)
likelihood = mash_compute_loglik(g,data)
weights=get_estimated_pi(g)


outname = head(labels,n=1)
write.table(c(likelihood),file = paste0('likelihoods/step',step,'/',outname,'.',cohort,'.',effect,'.txt'),row.names = FALSE,col.names = FALSE)
write.table(weights,file = paste0('likelihoods/step',step,'/',outname,'.',cohort,'.',effect,'.weights.txt'),col.names=FALSE)
