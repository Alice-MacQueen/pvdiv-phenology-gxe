## code to prepare `Taxa` dataset goes here

library(readr)
library(magrittr)
library(dplyr)

SNP_ID <- read_delim(file = file.path("data-raw",
                                      "Pvirgatum_4x_784g_imp_phased_maf0.019.012.indv"),
                     delim = "\n", col_names = FALSE)
Taxa <- SNP_ID %>%
  rename(PLANT_ID = X1)

save(Taxa, file = "data/Taxa.rda")
