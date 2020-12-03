library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)

#  devtools::install_github("alice-macqueen/switchgrassgwas")
inputfiles <- read_delim("~/Github/pvdiv-fitness-2019/analysis/inputkinship.txt", 
                         delim = " ", col_names = "SNPfiles")
outputdir <- file.path("~", "Github", "pvdiv-phenology-gxe", "analysis", "gwas", 
                       "Gulf_Midwest")

snp <- snp_attach(inputfiles$SNPfiles[4])
gwas_rds <- pvdiv_results_in_folder(path = outputdir, 
                                    pattern = "GWAS_datatable_(PKLE|STIL|KING|KBSM|CLMB|BRKG|LINC)_D2F") 
# Not enough Gulf & Midwest individuals at BRKG to run the univariate GWAS (<120).
phenotype_names <- str_sub(gwas_rds, start = 16, end = 23)
numSNPs <- ceiling(ceiling(1000000/length(phenotype_names)^2)/500)*500 # round 
# up number of SNPs to nearest 500.

mashin <- pvdiv_bigsnp2mashr(path = outputdir, 
                             snp = snp, gwas_rds = gwas_rds, numSNPs = numSNPs, 
                             phenotypes = phenotype_names, model = "linear", 
                             saveoutput = TRUE, scale = FALSE)

m <- mash_standard_run(path = outputdir, list_input = mashin,
                       saveoutput = TRUE, ref = NA)

source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
mash_plot_manhattan_by_condition(m = m, saveoutput = TRUE)
mash_plot_sig_by_condition(m = m, saveoutput = TRUE)


library(AnnotationDbi)
txdb <- loadDb(file = file.path("~", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))
annos <- pvdiv_table_topsnps(df = m, type = "mash", txdb = txdb, 
                             snp = snp)
saveRDS(annos, file = file.path(outputdir, paste0("Annotation_tables_mash_", 
                                                  "fourway_sites_", 
                                                  "Midwest_Gulf_D2F.rds")))

library(ashr)
for(n in 1:10){
  i <- get_significant_results(m)[n]
  
  effectplot <- switchgrassGWAS:::get_colnames(m) %>%
    enframe(name = NULL, value = "Conditions") %>%
    mutate(Conditions = str_sub(.data$Conditions, start = 6),
           mn = get_pm(m)[i,],
           se = get_psd(m)[i,])
  
  ggobject <- ggplot(data = effectplot) +
    geom_point(mapping = aes(x = as.factor(.data$Conditions), y = .data$mn)) +
    geom_errorbar(mapping = aes(ymin = .data$mn - .data$se,
                                ymax = .data$mn + .data$se,
                                x = .data$Conditions), width = 0.3) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Conditions", y = "Effect Size") +
    scale_x_discrete(labels = as_vector(.data$Conditions)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(filename = file.path(outputdir, paste0("Effect_plot_",
                              names(get_significant_results(m))[n], ".png")),
            plot = ggobject, base_aspect_ratio = 0.8, base_height = 4.5)
}

top_set <- mashin$top_set %>% mutate(start = POS, end = POS) %>% ungroup()

anno_tables <- pvdiv_table_topsnps(df = top_set, type = "table", 
                                   markers = markers, txdb = txdb, 
                                   anno_info = switchgrassGWAS::anno_info)

saveRDS(anno_tables[[1]], file.path(outputdir, paste0("Annotation_tables_", 
                                                      "mash_input_", numSNPs, 
                                                      "_SNPs_", 
                                                      "fourway_sites_", 
                                                      "Midwest_Gulf_D2F.rds")))
# ----------------------------
library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)

#  devtools::install_github("alice-macqueen/switchgrassgwas")
inputfiles <- read_delim("~/Github/pvdiv-fitness-2019/analysis/inputfull.txt", 
                         delim = " ", col_names = "SNPfiles")
outputdir <- file.path("~", "Github", "pvdiv-fitness-2019", "analysis", "gwas", 
                       "Full_log_biomass_survival")

snp <- snp_attach(inputfiles$SNPfiles[1])
gwas_rds <- pvdiv_results_in_folder(path = outputdir, 
                                    pattern = "GWAS_datatable_*.*BIOMASS") 

phenotype_names <- str_sub(gwas_rds, start = 16, end = -43)[1:11]
numSNPs <- ceiling(ceiling(1000000/length(phenotype_names)^2)/500)*500 # round 
# up number of SNPs to nearest 500.


mashin <- pvdiv_bigsnp2mashr(path = outputdir, 
                             snp = snp, gwas_rds = gwas_rds, numSNPs = numSNPs, 
                             phenotypes = phenotype_names, model = "linear", 
                             saveoutput = TRUE, scale = FALSE)

mash_output <- mash_standard_run(path = outputdir, list_input = mashin, 
                                 saveoutput = TRUE, ref = NA)

source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
mash_plot_manhattan_by_condition(m = mash_output, saveoutput = TRUE)
mash_plot_sig_by_condition(m = mash_output, saveoutput = TRUE)
mash_plot_effects(m = mash_output, n = 1, saveoutput = TRUE)

library(AnnotationDbi)
txdb <- loadDb(file = file.path("~", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))
annos <- pvdiv_table_topsnps(df = mash_output, type = "mash", txdb = txdb, 
                             snp = snp)
# ----------------------------


# -----------------------------------------------------------------
# Overlapping annotations for top SNPs, using raw log10p values
library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
library(AnnotationDbi)
txdb <- loadDb(file = file.path("~", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))
inputfiles <- read_delim("~/Github/pvdiv-fitness-2019/analysis/inputfull.txt", 
                         delim = " ", col_names = "SNPfiles")

snp <- snp_attach(inputfiles$SNPfiles[1])
markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos)

outputdir <- file.path("~", "Github", "pvdiv-fitness-2019", "analysis", "gwas", 
                       "Full_log_biomass_survival")

gwas_rds <- pvdiv_results_in_folder(path = outputdir, 
                                    pattern = "GWAS_datatable_") 

phenotype_names <- str_sub(gwas_rds, start = 16, end = -43)

mashin <- pvdiv_bigsnp2mashr(path = outputdir, snp = snp, gwas_rds = gwas_rds, 
                             numSNPs = 1000, phenotypes = phenotype_names, 
                             model = "linear", saveoutput = TRUE, clump = FALSE)
top_set <- mashin$top_set %>% mutate(start = POS, end = POS) %>% ungroup()

anno_tables <- pvdiv_table_topsnps(df = top_set, type = "table", 
                                   markers = markers, txdb = txdb, 
                                   anno_info = switchgrassGWAS::anno_info)

saveRDS(anno_tables[[1]], file.path(outputdir, paste0("Annotation_table_",
                                                 "full_set_log_biomass_survival_",
                                                 "10_sites_top1000_per_site", 
                                                 "linked_SNPs.rds")))