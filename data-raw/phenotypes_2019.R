# code to generate phenotypes_wide_2019 goes here

library(tidyverse)
library(XLConnect)
library(rlang)
library(switchgrassGWAS)
library(cowplot)
library(viridis)
library(bigsnpr)
library(scales)
library(lubridate)

source("../../Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
theme_set(theme_poster)

# Load all data

load("data/sites.rda")
load("data/metadata.rda")
load("data/Taxa.rda")

usa <- map_data("state")
wb_pheno <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                   paste0("GWAS DOE Switchgrass, ",
                                          "Consolidated Data"),
                                   paste0("DOE_GWAS_2019_Consolidated ",
                                          "Pheno Data_12-25-19.xlsx")))
lst_pheno <- readWorksheet(wb_pheno, sheet = getSheets(wb_pheno))
phenos_planting1 <- lst_pheno$`GWAS 2019 Greenup Data`
phenos_planting2 <- lst_pheno$`GWAS 2019 Core Phenotype Data`
phenos_key <- lst_pheno$`Column Key and Notes`

wb_KING2019 <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                      paste0("GWAS DOE Switchgrass, ",
                                             "Consolidated Data"),
                                      "Individual Site Files",
                                      "KING_GWAS_2019_Chlorosis Data.xlsx"))
lst_KING2019 <- readWorksheet(wb_KING2019, sheet = getSheets(wb_KING2019))
KING2019 <- as_tibble(lst_KING2019$`KING GWAS 2019 Chlorosis`)
key_KING19 <- as_tibble(lst_KING2019$`COLUMN KEY`)

wb_DEAD2018 <- loadWorkbook(file.path("C:", "Users", "ahm543", "Dropbox",
                                      paste0("GWAS DOE Switchgrass, ",
                                             "Consolidated Data"),
                                      "Forms and Miscellaneous",
                                      paste0("GWAS_2019_Mortality, Weak and ",
                                             "Not P. Virgatum Plants_ARCHIVE.xlsx")))
lst_DEAD2018 <- readWorksheet(wb_DEAD2018, sheet = getSheets(wb_DEAD2018))
DEAD2018 <- as_tibble(lst_DEAD2018$`GWAS 2019 Dead Only`)
LOST2018 <- as_tibble(lst_DEAD2018$`GWAS 2019 Dead, Weak, Not PV`)

LOST2018 <- LOST2018 %>%
  mutate(DEAD_2018 = case_when(
    NOTES %in% c("dead", "Dead", "Dead?", "DEAD") ~ 1,
    NOTES %in% c("1 tiller", "2 tillers", "3 tillers", "4 tillers",
                 "5 tillers", "Very Weak", "small crown", "Partial crown death",
                 "small crown/2 tillers", "small crown/3 tillers",
                 "small crown/4 tillers", "small crown/5 tillers",
                 "small crown/7 tillers", "small crown/8 tillers") ~ 0,
    NOTES %in% c("Not Virgatum", "switchgrass?", "Not Virgatum, Replant") ~ 2,
    TRUE ~ 0
  ))


phenos_2019 <- phenos_planting1 %>%
  left_join(LOST2018, by = c("SITE", "PLOT_GL", "PLOT_LC", "PLANT_ID",
                             "PLANT_.ID_GL")) %>%
  left_join(phenos_planting2, by = c("SITE", "PLOT_GL", "PLOT_LC")) %>%
  left_join(KING2019, by = "PLOT_GL") %>%
  mutate(PLANT_ID = ifelse(ACC == "AP13",
                           "AP13",
                           PLANT_ID)) %>%
  dplyr::select(-TC_EOS_2018, -`PLANT_.ID_GL`, -INDV, -(LD_EOS:BIOMASS))

# Drop 31 plants with confusion on whether plant died for now.

phenotypes_2019 <- phenos_2019 %>%
  filter(!(DEAD_2018 %in% c(1,2)) | is.na(DEAD_2018)) %>% # 1 = dead in 2018, and 2 = not a virgatum plant. Drop 1's and 2's. Keep 0's and NA's.
  filter(SRV != "N") %>% # Drop individuals that did not survive
  dplyr::select(-VERIFY, -PLOT_GL_CHK)

## --------------------------------
# Ensure all of the columns have the correct data type - mostly, numeric.
# Make time to event columns

phenotypes_2019 <- phenotypes_2019 %>%
  mutate(FRZ_DMG = as.numeric(FRZ_DMG),
         GR1 = as.numeric(GR1),
         GR50 = as.numeric(GR50),
         GR100 = as.numeric(GR100),
         EM1 = as.numeric(EM1),
         EM50 = as.numeric(EM50),
         FL1 = as.numeric(FL1),
         FL50 = as.numeric(FL50),
         TC_EOS = as.numeric(TC_EOS),
         HT_PAN_EOS = as.numeric(HT_PAN_EOS),
         DEAD_2018 = ifelse(is.na(DEAD_2018), 0, DEAD_2018),
         CHLR_101 = as.numeric(CHLR_101),
         CHLR_121 = as.numeric(CHLR_121),
         CHLR_136 = as.numeric(CHLR_136),
         CHLR_JASON_138 = as.numeric(CHLR_JASON_138),
         CHLR_150. = as.numeric(CHLR_150.),
         SPAD_158 = as.numeric(SPAD_158),
         SPADCOR_158 = as.numeric(SPADCOR_158),
         T_GR_EM = EM50 - GR50,
         T_GR_FL = FL50 - GR50,
         T_EM_FL = FL50 - EM50,
         EM50 = ifelse(EM50 - EM1 < 0, NA, EM50),
         FL50 = ifelse(FL50 - FL1 < 0, NA, FL50)
  )

# Get the mean phenotype value for each PLANT_ID at each SITE as the GWAS value
phemean_2019 <- phenotypes_2019 %>%
  group_by(PLANT_ID, SITE) %>%
  summarise(FRZ_DMG = mean(FRZ_DMG, na.rm = TRUE),
            GR1 = mean(GR1, na.rm = TRUE),
            GR50 = mean(GR50, na.rm = TRUE),
            GR100 = mean(GR100, na.rm = TRUE),
            EM1 = mean(EM1, na.rm = TRUE),
            EM50 = mean(EM50, na.rm = TRUE),
            FL1 = mean(FL1, na.rm = TRUE),
            FL50 = mean(FL50, na.rm = TRUE),
            TC_EOS = mean(TC_EOS, na.rm = TRUE),
            HT_PAN_EOS = mean(HT_PAN_EOS, na.rm = TRUE),
            T_GR_EM = mean(T_GR_EM, na.rm = TRUE),
            T_GR_FL = mean(T_GR_FL, na.rm = TRUE),
            T_EM_FL = mean(T_EM_FL, na.rm = TRUE),
            CHLR_101 = mean(CHLR_101, na.rm = TRUE),
            CHLR_121 = mean(CHLR_121, na.rm = TRUE),
            CHLR_136 = mean(CHLR_136, na.rm = TRUE),
            CHLR_JASON_138 = mean(CHLR_JASON_138, na.rm = TRUE),
            CHLR_150. = mean(CHLR_150., na.rm = TRUE),
            SPAD_158 = mean(SPAD_158, na.rm = TRUE),
            SPADCOR_158 = mean(SPADCOR_158, na.rm = TRUE)
  ) %>%
  mutate_all( ~ case_when(!is.nan(.x) ~ .x)) %>% # if the column is not NA, keep the value; else replace with NA as nothing is provided as another choice in case_when.
  filter(!is.na(SITE) & !is.na(PLANT_ID)) %>%
  left_join(metadata, by = "PLANT_ID") %>%
  dplyr::select(PLANT_ID:SITE, LIB_PHENO, Ecotype, GR1:SPADCOR_158, JL_genetic_subpop) %>%
  filter(LIB_PHENO == "Y")

phenotypes_wide_2019 <- phemean_2019 %>%
  dplyr::select(-LIB_PHENO, -Ecotype, -JL_genetic_subpop) %>%
  gather(c(GR1:SPADCOR_158), key = "Phenotype", value = value) %>%
  filter(!is.na(value)) %>%
  # saveRDS(file = "Cleaned pvdiv Phenotypes from 2018 for 784g 2019-05-16.rds")
  unite(col = "SITE_PHE", SITE, Phenotype) %>%
  spread(key = SITE_PHE, value = value) %>%
  right_join(Taxa) #%>% filter(PLANT_ID == "AP13")

save(phenotypes_wide_2019, file = "data/phenotypes_2019.rda")
