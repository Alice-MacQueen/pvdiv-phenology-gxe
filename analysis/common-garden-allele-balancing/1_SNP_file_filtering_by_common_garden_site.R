library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)


pvdiv_pop_subset_bigsnp <- function(subset, selector){
  subset_name <- paste0("Pvirgatum_4x_784g_imp_phased_", subset, "_subset")
  snpfile_subset <- subset(snp, ind.row = selector)
  snp_subset <- snp_attach(snpfile_subset)
  snp_writeBed(snp_subset, bedfile = paste0(subset_name, ".bed"))
  mac <- round(25/length(selector), digits = 2) # what should the MAF be?
  if(mac > 0.05){
    mac <- 0.05
  }
  plink <- download_plink()
  bedfile <- snp_plinkQC(plink, prefix.in = subset_name,
                         maf = mac, geno = 0.1, hwe = 1e-50, 
                         prefix.out = paste0(subset_name, "_maf", mac),
                         extra.options = "--allow-no-sex --allow-extra-chr --chr Chr01K Chr01N Chr02K Chr02N Chr03K Chr03N Chr04K Chr04N Chr05K Chr05N Chr06K Chr06N Chr07K Chr07N Chr08K Chr08N Chr09K Chr09N")
  rdsfile <- snp_readBed(bedfile)
}


# Load data - if you get errors, make sure to run R from the GWAS/ folder in pubjuenger
snp <- snp_attach("/home/pubjuenger/GWAS/Pvirgatum_4x_784g_imp_phased_maf0.005.rds")
switchgrassGWAS::phenotypes

# Use phenotypes to set up selection criteria for subsets. These need to be vectors of the row numbers in the Taxa dataframe that you have phenotypes for.
gwasct_sel <- which(!is.na(phenotypes$GWAS_CT)) # gwas_ct GWAS

## Run function to generate subset SNP files for the subset
pvdiv_pop_subset_bigsnp(subset = "gwasct", selector = gwasct_sel)
