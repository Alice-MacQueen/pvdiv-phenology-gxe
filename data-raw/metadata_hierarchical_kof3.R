library(tidyverse)

metadata732 <- read_csv("~/Github/pvdiv-climate-gwas/data/Metadata_2020-02.csv")
load("~/Github/pvdiv-phenology-gxe/data/metadata.rda")
ind_needed <- read_csv("~/Github/pvdiv-climate-gwas/results-33M-SNPs/Individuals_needed_for_climate_and_phenotype_GWAS.csv")

metadata732 <- metadata732 %>%
  add_row(PLANT_ID = "AP13", LIBRARY = "REF", Genetic_subpop_50per = "Texas",
          Genetic_subpop_95per = "Texas", Texas_Q = 1, Eastcoast_Q = 0,
          Midwest_DA = 0, LIB_GROWN = "Y", include_envGWAS = "N") %>%
  dplyr::select(PLANT_ID:LIB_GROWN)

metadata <- metadata %>%
  left_join(ind_needed) %>%
  dplyr::select(PLANT_ID, LIBRARY, LIB_GROWN, everything()) %>%
  mutate(LIB_GROWN = ifelse(is.na(LIB_GROWN), "N", LIB_GROWN),
         LIB_BIOCLIM = ifelse(is.na(include_clim), "N", include_clim)) %>%
  dplyr::select(PLANT_ID, LIBRARY, LIB_GROWN, LIB_BIOCLIM, everything())

metadata %>%
  left_join(metadata732) %>%
  unique() %>%
  dplyr::select(PLANT_ID, LIBRARY, LIB_GROWN, LIB_BIOCLIM, Genetic_subpop_50per:Midwest_DA, Ecotype, Ecotype_PC1, UP_ecotype:CO_ecotype, everything()) %>%
  dplyr::select(-(JL_genetic_subpop:include_clim)) %>%
  write_csv(., path = "data/metadata_hierarchical_kof3.csv")
