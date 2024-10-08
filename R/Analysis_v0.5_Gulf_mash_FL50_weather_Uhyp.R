# Gulf mash on FL50 non-interactive

## Load data
library(tidyverse)
library(bigsnpr)
library(mashr)
library(switchgrassGWAS)
library(rlist)
source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")

workingdir <- file.path("~", "Github", "pvdiv-phenology-gxe")
outputdir <- file.path(workingdir, "analysis", "gwas", "BLUPs", "single_sites")
metadata <- read_rds(file.path(workingdir, "data", "metadata.rds"))
sites <- read_rds(file.path(workingdir, "data", "sites.rds"))
subpops <- read_rds(file.path(workingdir, "data", "subpop_color_coding.rds"))
phe_gr_df <- read_rds(file.path(workingdir, "data",
                                "Weather_related_phe_for_asreml_h2_GR50.rds"))
phe_fl_df <- read_rds(file.path(workingdir, "data",
                                "Weather_related_phe_for_asreml_h2_FL50.rds"))
k_full <- read_rds("~/Github/pvdiv-genome/tensite_twoyear/Kinship_van_Raden_630_individuals_SNPs_r2_20percent.rds")

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

## Make U_hyp
U_hyp <- list()
ggcov <- list()
phe_hyp <- c("FL50", "dyln_fl50", "dyln_change_sec", "cgdd_12c_gr2fl", "crain_gr2fl", "crain_1d", "crain_2d", "crain_3d", "crain_5d", "crain_7d")

for(i in seq_along(phe_hyp)){
  for(k in seq_along(subpop_v2)){

    cov_df <- phe_fl_df %>%
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

    #ggcov1 <- pvdiv_plot_cov(cor_phe) + ggtitle(label = paste0(phe_hyp[i], "_", names(subpop_v2)[k]))
    #ggcov <- list.prepend(ggcov, ggcov1)
    U_hyp <- list.prepend(U_hyp, cor_phe)
    names(U_hyp)[1] <- paste0(phe_hyp[i], "_", names(subpop_v2)[k])
  }
}

# --------------------
## Run mash
for(k in c(3)){
  phe_suffix <- paste("FL50", subpop_v[k],
                      sep = "*.*")
  gwas_rds <- pvdiv_results_in_folder(path = file.path(outputdir),
                                      pattern = paste0("GWAS_datatable*.*",
                                                       phe_suffix))
  phenotypes_1 <- str_sub(gwas_rds, start = 16, end = -42) %>%
    str_replace(., "_$", "")
  phenotypes <- paste0(names(subpop_v)[k], "_", phenotypes_1)
  numSNPs <- round(ceiling(1200000/length(phenotypes)^2)/1000)*1000
  U_ed <- readRDS(file.path(outputdir, paste0("Mash_data_driven_covariances_",
                                              numSNPs, "_FL50_",
                                              names(subpop_v)[k], ".rds")))
  m_out <- mash_standard_run(path = outputdir, numSNPs = numSNPs,
                             U_hyp = U_hyp, U_ed = U_ed,
                             suffix = paste0("FL50_", names(subpop_v)[k]),
                             saveoutput = TRUE)
}
