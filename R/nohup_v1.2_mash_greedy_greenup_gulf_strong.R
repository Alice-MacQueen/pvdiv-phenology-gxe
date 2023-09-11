library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
library(mashr)
library(rlist)
library(here)

source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
outputdir <- here("analysis", "mash", "Gulf", "inputs")
mashdir <- here("analysis", "mash_greedy_algorithm")

#metadata <- read_rds(here("data", "metadata.rds"))
# pvdiv PLANT_ID metadata
#sites <- read_rds(here("data", "sites.rds"))
# pvdiv common garden metadata
#subpops <- read_rds(here("data", "subpop_color_coding.rds"))
phe_gr_df <- read_rds(here("data",
                           "Weather_related_phe_for_asreml_h2_GR50.rds"))
# all greenup related phenotypes, including greenup as functions of weather variables.
phe_fl_df <- read_rds(here("data",
                           "Weather_related_phe_for_asreml_h2_FL50.rds"))
# all flowering related phenotypes, including functions of weather variables.
k_full <- read_rds(here("data",
                        "Kinship_van_Raden_630_individuals_SNPs_r2_20percent.rds"))
# kinship matrix for 630 tetraploid individuals.

## vectors with subsets of gardens, phenotypes, and genetic subpopulations.
#site_v <- list(TX1 = "TX1", TX2 = "TX2", TX3 = "TX3", OK = "OK", MO = "MO",
#               NE = "NE", MI = "MI", SD = "SD")#, eight_sites = c("TX1", "TX2", "TX3", "OK", "MO", "NE", "MI", "SD"), tx_sites = c("TX1", "TX2", "TX3"), north_sites = c("MO", "NE", "MI", "SD"))
subpop_v <- list(Gulf_and_Midwest = "363g_12.1M", Midwest = "134g_8.8M", 
                 Gulf = "229g_10.3M")
#subpop_v2 <- list(Gulf_and_Midwest = c("Gulf", "Midwest"), Midwest = "Midwest", 
#                  Gulf = "Gulf")
#phe_single <- c("FL50", "dyln_fl50", "cgdd_12c_gr2fl", "crain_gr2fl",
#                "crain_1d", "GR50", "cgdd_12c_10d")
#phe_mult <- list(flowering_5 = c("FL50", "dyln_fl50", "cgdd_12c_gr2fl",
#                                 "crain_gr2fl", "crain_1d"), 
#                 greenup_2 = c("GR50", "cgdd_12c_10d"))
#inputfiles <- read_delim(here("analysis", "gwas",
#                              "inputkinship.txt"),  delim = " ", 
#                         col_names = "SNPfiles") 
subpop_v3 <- list(Gulf_and_Midwest = "all", Midwest = "midwest", Gulf = "gulf")
# use this one for the likelihood paths
# Read in likelihood paths, then use that to read in matrices, then use those to run mash.

stopcondpath <- here("analysis", "mash_greedy_algorithm", "mash_switchgrass_models", "mash", 
                     "all_matrices", "likelihood_paths")
U_greedy <- read_csv(file = file.path(stopcondpath, 
                                      paste0("greenup", ".", "gulf", 
                                             ".strong.stop.condition.csv"))) %>%
  filter(matrix != "0.0")
## Read in Uhyp
Upath <- here("analysis", "mash_greedy_algorithm", "mash_switchgrass_models", 
              "mash", "all_matrices", "matrices")
U_list <- list()
for (i in 1:nrow(U_greedy)) {
  U_single <- read_csv(file = file.path(Upath, U_greedy$matrix[i]))
  U_single <- U_single %>%
    select(-`...1`)
  U_single <- as.matrix(unname(U_single))
  U_list <- list.append(U_list, U_single)
}
names(U_list) <- str_remove(U_greedy$matrix, ".csv")


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

list_input <- switchgrassGWAS:::load_mash_df(path = outputdir, numSNPs = numSNPs,
                                             suffix = suffix)
# Just run mash itself as the greedy choice of cov matrices changes mash run
# enough that mash_standard_run will not perform correctly (eg it always adds
# all canonical covariance matrices).
Bhat_strong <- as.matrix(list_input$B_hat_strong)
Shat_strong <- as.matrix(list_input$S_hat_strong)
Bhat_random <- as.matrix(list_input$B_hat_random)
Shat_random <- as.matrix(list_input$S_hat_random)


data_r <- mashr::mash_set_data(Bhat_random, Shat_random)
Vhat <- mashr::estimate_null_correlation_simple(data = data_r)

data_strong <- mashr::mash_set_data(Bhat_strong, Shat_strong, V = Vhat)
data_random <- mashr::mash_set_data(Bhat_random, Shat_random, V = Vhat)

U_single <- as.matrix(U_single)
# Run mash on the random dataset using the random data w/ correlation structure
message(paste0("Fit mash to the random tests using both data-driven and ",
               "canonical covariances."))
U_c <- mashr::cov_canonical(data_random)
m = mashr::mash(data_random, Ulist = U_list, outputlevel = 1)

# Run mash on the strong dataset (or all data) using
# the previous results from the random data
message(paste0("Compute posterior matrices for the strong effects",
               " using the mash fit from the
                 random tests."))
m2 = mashr::mash(data_strong, g = get_fitted_g(m), fixg = TRUE)

saveRDS(m2, file = file.path(mashdir, names(subpop_v)[k], 
                                 paste0("Strong_Effects_", numSNPs, "_SNPs_",
                                        suffix, "_Uhyp.rds")))

print("Log likelihood with specified covariance matrices: ")
print(get_loglik(m2), digits = 10)
print("How many significant markers?")
print(length(get_significant_results(m2)))
switchgrassGWAS::mash_plot_covar(m2)
switchgrassGWAS::mash_plot_manhattan_by_condition(m2)
switchgrassGWAS::mash_plot_pairwise_sharing(m2)
gxe <- switchgrassGWAS::get_GxE(m2)
