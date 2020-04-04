library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
library(AnnotationDbi)
library(multtest)
library(viridis)
library(cowplot)
library(data.table)
source("~/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
get_date_filename <- function(){
  str_replace_all(str_replace_all(Sys.time(), ":", "."), " ", "_")
}


# ---------- Choose SNP file & SVD/Kinship matrix appropriate for number of sites
snp <- snp_attach(file.path("~", "Github", "pvdiv-genome",
                            "tensite_twoyear",
                            "Pvirgatum_V5_GWAS_630g_33M_tensite_twoyear.rds"))
svd <- readRDS(file.path("~", "Github", "pvdiv-genome",
                         "tensite_twoyear",
                         "SVD_20PCs_630g_33M_tensite_twoyear.rds"))

# ---------- Load or make some useful objects ------

NCORES <- nb_cores()
G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
plants <- snp$fam$sample.ID
bonferroni <- -log10(0.05/length(POS))
markers <- tibble(CHR = CHR, POS = POS)

txdb <- loadDb(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

anno_info <-
  read_delim(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                              "Pvirgatum_516_v5.1.annotation_info.txt"),
             col_names = TRUE, delim = "\t")
gene_anno_info <- tbl_df(anno_info) %>%
  distinct(locusName, .keep_all = TRUE)


# ------- Set analysis directory and phenotypes
analysisdir <- file.path("~", "Github", "pvdiv-phenology-gxe", "analysis",
                         "gwas", "linear")
setwd(analysisdir)
phe_df <- read_rds(file.path("~", "Github",
                             "pvdiv-phenotypes", "data",
                             "Phenotypes_D2F_fourway_project.rds"))


df <- plants %>%
  enframe(name= NULL, value = "PLANT_ID") %>%
  left_join(phe_df, by = "PLANT_ID")

for(i in 2:nrow(df)){

  # ------- Choose best pop structure correction --------------------------

  df1 <- df %>%
    dplyr::select(PLANT_ID, all_of(i))

  lambdagc <- pvdiv_lambda_GC(df = df1, type = "linear", snp = snp,
                              covar = svd, ncores = NCORES,
                              npcs = c(0:20), saveoutput = TRUE)


  # ------- Run GWAS with best pop structure correction -------------------
  # Use a new function asr_best_PC_df(df) to find the NumPCs where lambda_GC
  # is closest to one.

  PCdf <- switchgrassGWAS:::pvdiv_best_PC_df(lambdagc) # asv_best_PC_df(lambdagc)
  PCdf1 <- PCdf[1,]
  # PCdf1 <- data.frame(NumPCs = 20)

  if(PCdf1$NumPCs == 0){
    gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, ncores = NCORES)
  } else {
    gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, covar = svd,
                       ncores = NCORES, npcs = PCdf1$NumPCs)
  }

  savegwas = TRUE
  # Save a data.table object with the GWAS results
  gwas_data <- data.table(CHR = markers$CHR, POS = markers$POS,
                          estim = gwas$estim, std_err = gwas$std.err,
                          bigsnpscore = gwas$score)
  gwas_data[,pvalue := predict(gwas, log10 = FALSE)]
  gwas_data[,log10p := -log10(pvalue)]
  gwas_data[,FDR_adj := p.adjust(pvalue, method = "BH")]
  if(savegwas == TRUE){
  write_rds(gwas_data, path = paste0("GWAS_datatable_", names(df1)[2], "_",
                                     PCdf1$NumPCs, "_PCs", "_.rds"),
            compress = "gz")
  }

  saveplots = TRUE
  if(saveplots == TRUE){
  # Find 10% FDR threshold
  FDRthreshhi <- gwas_data %>%
    as_tibble() %>%
    filter(between(FDR_adj, 0.10001, 1)) %>%  # 10% FDR threshold currently.
    summarise(thresh = max(log10p))
  FDRthreshlo <- gwas_data %>%
    as_tibble() %>%
    filter(between(FDR_adj, 0, 0.09999)) %>%  # 10% FDR threshold currently.
    summarise(thresh = min(log10p))
  if(FDRthreshhi$thresh[1] > 0 & FDRthreshlo$thresh[1] > 0){
    FDRthreshold = (FDRthreshhi$thresh[1] + FDRthreshlo$thresh[1])/2
  } else if(FDRthreshhi$thresh[1] > 0){
    FDRthreshold = FDRthreshhi$thresh[1]
  } else if(FDRthreshlo$thresh[1] > 0){
    FDRthreshold = FDRthreshlo$thresh[1]
  } else {
    FDRthreshold = NA
  }
  # Save a Manhattan plot with 10% FDR
  ggmanobject1 <- gwas_data %>%
    filter(log10p > 1) %>%
    ggplot(aes(x = POS, y = log10p)) +
    geom_hline(yintercept = c(5, 10), color = "lightgrey") +
    geom_point(aes(color = CHR, fill = CHR)) +
    geom_hline(yintercept = FDRthreshold, color = "black", linetype = 2,
               size = 1) +
    facet_wrap(~ CHR, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_color_viridis(option = "B", end = 0.8, discrete = TRUE) +
    scale_fill_viridis(option = "B", end = 0.8, discrete = TRUE) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA), legend.position = "none") +
    labs(x = "Chromosome", y = "-log10(p value)") +
    scale_x_continuous(expand = c(0.18, 0.18))

  save_plot(paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs, "_PCs_",
                   "10percent_FDR_", get_date_filename(),
                   ".png"), plot = ggmanobject1, base_asp = 4, base_height = 4)

  # Save a QQplot
  ggqqplot <- pvdiv_qqplot(ps = gwas_data$pvalue, lambdaGC = TRUE)
  save_plot(paste0("QQplot_", names(df1)[2], "_", PCdf1$NumPCs, "_PCs_FDR_",
                   get_date_filename(), ".png"),
            plot = ggqqplot, base_asp = 1, base_height = 4)

  # Save a Manhattan plot with Bonferroni
  ggmanobject2 <- gwas_data %>%
    filter(log10p > 1) %>%
    ggplot(aes(x = POS, y = log10p)) +
    geom_hline(yintercept = c(5, 10), color = "lightgrey") +
    geom_point(aes(color = CHR, fill = CHR)) +
    geom_hline(yintercept = bonferroni, color = "black", linetype = 2,
               size = 1) +
    facet_wrap(~ CHR, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_color_viridis(option = "B", end = 0.8, discrete = TRUE) +
    scale_fill_viridis(option = "B", end = 0.8, discrete = TRUE) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA), legend.position = "none") +
    labs(x = "Chromosome", y = "-log10(p value)") +
    scale_x_continuous(expand = c(0.18, 0.18))

  save_plot(paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs,
                   "_PCs_Bonferroni_", get_date_filename(),
                   ".png"), plot = ggmanobject2, base_asp = 4, base_height = 4)
  }

  saveannos = TRUE
  if(saveannos == TRUE){
  ## Save annotation tables for the top associations
  anno_tables <- pvdiv_table_topsnps(df = gwas, type = "bigsnp", n = c(10,100),
                                     FDRalpha = 0.1,
                                     rangevector = c(0, 50000),
                                     markers = markers,
                                     anno_info = gene_anno_info,
                                     txdb = txdb)
  saveRDS(anno_tables, file.path("/", "home", "alice", "Github",
                                 "pvdiv-phenology-gxe", "analysis",
                                 "all-phenotypes-two-years",
                                 paste0("Annotation_tables_", names(df1)[2],
                                        "_", PCdf1$NumPCs, "_PCs", ".rds")))
  }

}
