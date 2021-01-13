library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
library(mashr)
library(rlist)

source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
workingdir <- file.path("~", "Github", "pvdiv-phenology-gxe")
outputdir <- file.path(workingdir, "analysis", "gwas", "BLUPs", "single_sites")
mashdir <- file.path(workingdir, "analysis", "mash")

metadata <- read_rds(file.path(workingdir, "data", "metadata.rds"))
# pvdiv PLANT_ID metadata
sites <- read_rds(file.path(workingdir, "data", "sites.rds"))
# pvdiv common garden metadata
subpops <- read_rds(file.path(workingdir, "data", "subpop_color_coding.rds"))
phe_gr_df <- read_rds(file.path(workingdir, "data",
                                "Weather_related_phe_for_asreml_h2_GR50.rds"))
# all greenup related phenotypes, including greenup as functions of weather variables.
phe_fl_df <- read_rds(file.path(workingdir, "data",
                                "Weather_related_phe_for_asreml_h2_FL50.rds"))
# all flowering related phenotypes, including functions of weather variables.
k_full <- read_rds(file.path(workingdir, "data",
                             "Kinship_van_Raden_630_individuals_SNPs_r2_20percent.rds"))
# kinship matrix for 630 tetraploid individuals.

## vectors with subsets of gardens, phenotypes, and genetic subpopulations.
site_v <- list(TX1 = "TX1", TX2 = "TX2", TX3 = "TX3", OK = "OK", MO = "MO",
               NE = "NE", MI = "MI", SD = "SD")#, eight_sites = c("TX1", "TX2", "TX3", "OK", "MO", "NE", "MI", "SD"), tx_sites = c("TX1", "TX2", "TX3"), north_sites = c("MO", "NE", "MI", "SD"))
subpop_v <- list(Gulf_and_Midwest = "363g_12.1M", Midwest = "134g_8.8M", 
                 Gulf = "229g_10.3M")
subpop_v2 <- list(Gulf_and_Midwest = c("Gulf", "Midwest"), Midwest = "Midwest", 
                  Gulf = "Gulf")
phe_single <- c("FL50", "dyln_fl50", "cgdd_12c_gr2fl", "crain_gr2fl",
                "crain_1d", "GR50", "cgdd_12c_10d")
phe_mult <- list(flowering_5 = c("FL50", "dyln_fl50", "cgdd_12c_gr2fl",
                                 "crain_gr2fl", "crain_1d"), 
                 greenup_2 = c("GR50", "cgdd_12c_10d"))
inputfiles <- read_delim(file.path(workingdir, "analysis", "gwas",
                                   "inputkinship.txt"),  delim = " ", 
                         col_names = "SNPfiles") 
# paths to large SNP files


## Define Uhyp

U_hyp <- list()
phe_hyp <- c("cgdd_12c_10d", "cgdd_12c_18d", "cgdd_12c_5d", 
             "tave_5d", "tave_10d", "tave_18d")
#phe_hyp <- c("dyln_fl50", "dyln_change_sec", "cgdd_12c_gr2fl", 
#             "crain_gr2fl", "crain_1d", "crain_3d", "crain_5d")

for(i in seq_along(phe_hyp)){
  for(k in seq_along(subpop_v2)){
    
    cov_df <- phe_gr_df %>%
      ungroup() %>%
      left_join(select(sites, SITE, manu_site), by = "SITE") %>%
      filter(SUBPOP %in% subpop_v2[[k]]) %>%
      select(PLANT_ID, SUBPOP, manu_site, phe_hyp[i]) %>%
      group_by(PLANT_ID, SUBPOP, manu_site) %>%
      mutate(IND = row_number()) %>%
      pivot_wider(names_from = manu_site, values_from = phe_hyp[i]) %>%
      mutate(PLANT_ID = as.factor(PLANT_ID),
             IND = as.factor(IND)) %>%
      select(-IL) %>%
      select(PLANT_ID, SUBPOP, IND, MI, MO, NE, OK, SD, TX1, TX2, TX3)
    
    cor_phe <- cor(cov_df[,-(1:3)], use = "pairwise")
    diag(cor_phe) <- cov_df %>%
      ungroup() %>%
      summarise(MI = sqrt(var(MI, na.rm = TRUE))/mean(MI, na.rm = TRUE),
                MO = sqrt(var(MO, na.rm = TRUE))/mean(MO, na.rm = TRUE),
                NE = sqrt(var(NE, na.rm = TRUE))/mean(NE, na.rm = TRUE),
                OK = sqrt(var(OK, na.rm = TRUE))/mean(OK, na.rm = TRUE),
                SD = sqrt(var(SD, na.rm = TRUE))/mean(SD, na.rm = TRUE),
                TX1 = sqrt(var(TX1, na.rm = TRUE))/mean(TX1, na.rm = TRUE),
                TX2 = sqrt(var(TX2, na.rm = TRUE))/mean(TX2, na.rm = TRUE),
                TX3 = sqrt(var(TX3, na.rm = TRUE))/mean(TX3, na.rm = TRUE)) %>%
      as_vector()

    U_hyp <- list.prepend(U_hyp, cor_phe)
    names(U_hyp)[1] <- paste0(phe_hyp[i], "_", names(subpop_v2)[k])
  }
}

# Files with the six suffixes I have/want mash results for:
numSNPs_v <- c(19000, 19000, 19000, 19000, 24000, 33000)
suffix_v <- c("FL50_Gulf_and_Midwest", "FL50_Gulf", "GR50_Gulf_and_Midwest",
              "GR50_Gulf", "GR50_Midwest", "FL50_Midwest")
i = 4  # change this for different scripts
# B_hat_strong_df_19000topSNPs_FL50_Gulf_and_Midwest.rds
# B_hat_strong_df_19000topSNPs_FL50_Gulf.rds
# B_hat_strong_df_19000topSNPs_GR50_Gulf_and_Midwest.rds
# B_hat_strong_df_19000topSNPs_GR50_Gulf.rds
# B_hat_strong_df_24000topSNPs_GR50_Midwest.rds
# B_hat_strong_df_33000topSNPs_FL50_Midwest.rds

## mash Gulf greenup
k = 3
numSNPs <- numSNPs_v[i]
suffix <- suffix_v[i]
U_ed <- file.path(mashdir, paste0("Mash_data_driven_covariances_", numSNPs, "_", 
                                  suffix, ".rds"))

m_out <- mash_standard_run(path = file.path(outputdir), 
                           U_ed = U_ed,
                           numSNPs = numSNPs, 
                           suffix = suffix,
                           saveoutput = FALSE)
saveRDS(m_out, file = file.path(mashdir, names(subpop_v)[k], 
                                paste0("Strong_Effects_", numSNPs, "_SNPs_",
                                       suffix, "_no_Uhyp.rds")))
m_out1 <- mash_standard_run(path = file.path(outputdir), 
                            U_ed = U_ed,
                            numSNPs = numSNPs, 
                            U_hyp = U_hyp,
                            suffix = suffix,
                            saveoutput = FALSE)
saveRDS(m_out1, file = file.path(mashdir, names(subpop_v)[k], 
                                 paste0("Strong_Effects_", numSNPs, "_SNPs_",
                                        suffix, "_Uhyp.rds")))
