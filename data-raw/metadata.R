## code to prepare `metadata` dataset goes here
library(tidyverse)
library(XLConnect)
wb_metadata <- loadWorkbook(file.path("data-raw",
                                      "PVDIV_Metadata_2019-09-03.xlsx"))
lst_metadata <- readWorksheet(wb_metadata, sheet = getSheets(wb_metadata))

pvdiv <- as_tibble(lst_metadata$`PVDIV MASTER METADATA`)

kof5 <- read_csv("data-raw/k3.nnet.gen.strict.ecotypeClassification.csv")

eco_geno_subgroups <- kof5 %>%
  dplyr::select(PLANT_ID:pc1_53.7perc, main.subpop:NEC) %>%
  dplyr::select(-train.id) %>%
  dplyr::rename(UP_ecotype = Upland, TX_ecotype = TX.lowland,
                CO_ecotype = Coast.lowland, Ecotype = eco.class,
                Ecotype_PC1 = pc1_53.7perc, JL_genetic_subpop = main.subpop,
                JL_detail_genetic = detail.subpop, TX_genetic = TX,
                GC_genetic = GC, MW_genetic = MW, SEC_genetic = SEC,
                NEC_genetic = NEC)

pvdiv <- pvdiv %>%
  select(PLANT_ID:LIB_NOTES, PLOIDY, COLLECTION_TYPE:COLLECTION_METHOD, COLLECTOR, STATE:COLL_DATE) %>%
  mutate(LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE),
         ELEVATION = as.numeric(ELEVATION),
         LIBRARY = ifelse(LIBRARY == "NA", NA, LIBRARY),
         LIB_CLIMATE = ifelse(LIB_CLIMATE == "NA", NA, LIB_CLIMATE),
         LIB_PHENO = ifelse(LIB_PHENO == "NA", NA, LIB_PHENO),
         LIB_ACC = ifelse(LIB_ACC == "NA", NA, LIB_ACC),
         COLL_DATE = ifelse(COLL_DATE == "NA", NA, COLL_DATE),
         HABITAT = ifelse(HABITAT == "NA", NA, HABITAT),
         LOCALITY = ifelse(LOCALITY == "NA", NA, LOCALITY),
         COUNTY = ifelse(COUNTY == "NA", NA, COUNTY),
         STATE = ifelse(STATE == "NA", NA, STATE),
         COLLECTOR = ifelse(COLLECTOR == "NA", NA, COLLECTOR),
         NOTE_LATLONG = ifelse(NOTE_LATLONG == "NA", NA, NOTE_LATLONG))

metadata <- pvdiv %>%
  full_join(eco_geno_subgroups, by = "PLANT_ID") %>%
  arrange(desc(LIB_PHENO), PLANT_ID)

save(metadata, file = "data/metadata.rda")
