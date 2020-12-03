library(bigsnpr)
library(tidyverse)
library(switchgrassGWAS)
library(AnnotationDbi)
txdb <- loadDb(file = file.path("~", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

outputdir <- file.path("~", "Github", "pvdiv-phenology-gxe", "analysis", "gwas", "GDD")

inputfiles <- read_delim("~/Github/pvdiv-phenology-gxe/analysis/gwas/inputkinship.txt", delim = " ", col_names = "SNPfiles") # Midwest & Gulf, Midwest, Gulf, Atlantic, Midwest & Atlantic, Gulf & Atlantic

phe_subsets <- read_rds("~/Github/pvdiv-phenology-gxe/data/Phenotype_Subsets_for_Subpop_GWAS.rds")

for(i in 1:nrow(inputfiles)){
  phe_df <- phe_subsets[[i]]
  snp <- snp_attach(inputfiles$SNPfiles[i])
  pvdiv_standard_gwas(snp, df = phe_df, type = "linear", outputdir = outputdir,
                      savegwas = TRUE, saveannos = FALSE, minphe = 100)
}
