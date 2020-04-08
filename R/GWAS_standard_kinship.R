library(bigsnpr)
library(tidyverse)
library(switchgrassGWAS)
library(AnnotationDbi)
txdb <- loadDb(file = file.path("~", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

outputdir <- file.path("~", "Github", "pvdiv-phenology-gxe", "analysis", "gwas")

inputfiles <- read_delim("~/Github/pvdiv-phenology-gxe/analysis/gwas/inputkinship.txt", delim = " ", col_names = "SNPfiles")

phe_df <- read_rds(file.path("~", "Github",
                             "pvdiv-phenology-gxe", "data",
                             "Phenotypes_D2F_fourway_project.rds"))
for(i in 2:nrow(inputfiles)){
  snp <- snp_attach(inputfiles$SNPfiles[i])
  pvdiv_standard_gwas(snp, df = phe_df, type = "linear", outputdir = outputdir,
                      savegwas = TRUE, saveannos = TRUE, txdb = txdb,
                      minphe = 120)
}
