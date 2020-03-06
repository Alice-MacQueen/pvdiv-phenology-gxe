library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")


NCORES <- nb_cores()

snp <- snp_attach(file.path("/", "home", "alice", "Github", "pvdiv-genome",
                            "tensite_twoyear",
                            "Pvirgatum_V5_GWAS_630g_33M_tensite_twoyear.rds"))

G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
plants <- snp$fam$sample.ID
bonferroni <- -log10(0.05/length(POS))

CHRN <- enframe(CHR, name = NULL) %>%
  dplyr::rename(CHR = .data$value) %>%
  mutate(CHRN = case_when(.data$CHR == "Chr01K" ~ 1,
                          .data$CHR == "Chr01N" ~ 2,
                          .data$CHR == "Chr02K" ~ 3,
                          .data$CHR == "Chr02N" ~ 4,
                          .data$CHR == "Chr03K" ~ 5,
                          .data$CHR == "Chr03N" ~ 6,
                          .data$CHR == "Chr04K" ~ 7,
                          .data$CHR == "Chr04N" ~ 8,
                          .data$CHR == "Chr05K" ~ 9,
                          .data$CHR == "Chr05N" ~ 10,
                          .data$CHR == "Chr06K" ~ 11,
                          .data$CHR == "Chr06N" ~ 12,
                          .data$CHR == "Chr07K" ~ 13,
                          .data$CHR == "Chr07N" ~ 14,
                          .data$CHR == "Chr08K" ~ 15,
                          .data$CHR == "Chr08N" ~ 16,
                          .data$CHR == "Chr09K" ~ 17,
                          .data$CHR == "Chr09N" ~ 18,
                          TRUE ~ 19
  ))
markers <- tibble(CHR = CHR, POS = POS)


# Run bigsnp2mashr

markers20per <- snp_clumping(G, infos.chr = CHRN$CHRN, ncores = NCORES)
markers2 <- markers[markers20per,] # Get unlinked markers for mash
gwas_rds <- pvdiv_results_in_folder(path = "~/Github/pvdiv-phenology-gxe/analysis/all-phenotypes-two-years/BIOMASS", pattern = "GWAS_object_") # 8 phenotypes from 4 sites - 32 phenotypes total with > 200 individuals phenotyped and sequenced.

phenotype_names <- str_sub(gwas_rds, start = 13, end = -5)

numSNPs <- ceiling(1000000/length(phenotype_names)^2)
numSNPs <- 12500 # the above gave 20409, so just rounded this to 20000.

mash_input <- pvdiv_bigsnp2mashr(path = "~/Github/pvdiv-phenology-gxe/analysis/all-phenotypes-two-years/BIOMASS", gwas_rds = gwas_rds,
                                 phenotypes = phenotype_names, numSNPs = numSNPs,
                                 markers = markers, markers2 = markers2, G = G,
                                 model = "linear", saveoutput = TRUE)

mash_output <- mash_standard_run(path = "~/Github/pvdiv-phenology-gxe/analysis/all-phenotypes-two-years/BIOMASS", list_input = mash_input,
                                 numSNPs = numSNPs, saveoutput = TRUE)

anno_tables <- pvdiv_table_topsnps(df = mash_output, type = "mash",
                                   n = c(10,500), FDRalpha = NA,
                                   rangevector = c(0, 20000, 100000),
                                   markers = markers,
                                   anno_info = gene_anno_info,
                                   txdb = txdb)
saveRDS(anno_tables, file.path("~", "Github", "pvdiv-phenology-gxe", "analysis", 
                               "all-phenotypes-two-years", "BIOMASS",
                               paste0("Annotation_tables_",
                                      "Flowering_date_50per_mash_",
                                      "7_sites_4_PCs.rds")))

mash_plot_manhattan_by_condition(m = m, saveoutput = TRUE)
mash_plot_sig_by_condition(m=m, saveoutput = TRUE)
mash_plot_effects(m=m, n=1, saveoutput = TRUE)