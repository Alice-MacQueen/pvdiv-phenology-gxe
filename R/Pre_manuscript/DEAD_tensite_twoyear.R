library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)
library(AnnotationDbi)
library(multtest)
library(viridis)
library(cowplot)
library(data.table)
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

svd <- readRDS(file.path("/", "home", "alice", "Github", "pvdiv-genome",
                         "tensite_twoyear", "SVD_20PCs_630g_33M_tensite_twoyear.rds"))
load("~/Github/pvdiv-phenology-gxe/data/phenotypes_2018.rda")

analysisdir <- "~/Github/pvdiv-phenology-gxe/analysis/DEAD_2018/"
setwd(analysisdir)
dead <- readRDS("Logical_SITE_DEAD_2018.rds")

txdb <- loadDb(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

anno_info <-
  read_delim(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                              "Pvirgatum_516_v5.1.annotation_info.txt"),
             col_names = TRUE, delim = "\t")
gene_anno_info <- tbl_df(anno_info) %>%
  distinct(locusName, .keep_all = TRUE)

df <- plants %>%
  enframe(value = "PLANT_ID") %>%
  left_join(phenotypes_2018, by = "PLANT_ID") %>%
  dplyr::select(-name, PLANT_ID,  ends_with("DEAD_2018"))
phenotypes_2018 %>%
  dplyr::select(PLANT_ID,  ends_with("DEAD_2018"))

PCdf <- read_csv("~/Github/pvdiv-phenology-gxe/analysis/all-phenotypes-two-years/PCdf2_BIOMASS.csv")

for(i in 1:nrow(PCdf)){
  PCdf1 <- PCdf[i,]
  df1 <- df %>%
    dplyr::select(PLANT_ID, PCdf1$trait)

  if(PCdf1$NumPCs == 0){
    gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, ncores = NCORES)
  } else {
    gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, covar = svd,
                       ncores = NCORES, npcs = PCdf1$NumPCs)
  }

  gwas_df <- gwas %>%
    mutate(p.value = predict(gwas, log10 = FALSE),
           log10p = -log10(p.value))
  res <- mt.rawp2adjp(gwas_df$p.value, alpha = 0.1, proc = "BH")
  adj_p <- res$adjp[order(res$index), ]
  gwas_df <- cbind(markers, gwas_df, adj_p) %>%
    as_tibble()

  gwas_data <- data.table(chr = markers$CHR, pos = markers$POS,
                          estim = gwas_df$estim, std_err = gwas_df$std.err,
                          bigsnpscore = gwas_df$score,
                          pvalue = gwas_df$p.value, log10p = gwas_df$log10p,
                          FDR_10per = gwas_df$BH)
  saveRDS(gwas_data, file = paste0("GWAS_datatable_", names(df1)[2], "_",
                                   PCdf1$NumPCs, "_PCs", "_.rds"))
  FDRthresh <- gwas_data %>%
    as_tibble() %>%
    filter(between(FDR_10per, 0.01, 0.1)) %>%
    summarise(thresh = max(log10p))

  ggmanobject1 <- gwas_data %>%
    filter(log10p > 1) %>%
    ggplot(aes(x = pos, y = log10p)) +
    geom_hline(yintercept = c(10,15), color = "lightgrey") +
    geom_point(aes(color = chr, fill = chr)) +
    geom_hline(yintercept = FDRthresh$thresh[1], color = "black", linetype = 2,
               size = 1) +
    facet_wrap(~ chr, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_color_viridis(option = "B", end = 0.8, discrete = TRUE) +
    scale_fill_viridis(option = "B", end = 0.8, discrete = TRUE) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA), legend.position = "none") +
    labs(x = "Chromosome", y = "-log10(p value)") +
    scale_x_continuous(expand = c(0.18, 0.18))
  #scale_shape_manual(values = rep(c(21,22),9), guide = FALSE)

  save_plot(paste0("Manhattan_", names(df1)[2], "_",
                   PCdf1$NumPCs, "_PCs_FDR_", str_replace_all(Sys.time(), ":", "."),
                   ".png"),
            plot = ggmanobject1, base_aspect_ratio = 4, base_height = 4)

  # Save annotation tables for the top associations
  anno_tables <- pvdiv_table_topsnps(df = gwas, type = "bigsnp", n = c(10,500),
                                     FDRalpha = NA,
                                     rangevector = c(0, 20000, 100000),
                                     markers = markers,
                                     anno_info = gene_anno_info,
                                     txdb = txdb)
  saveRDS(anno_tables, file.path("/", "home", "alice", "Github",
                                 "pvdiv-phenology-gxe", "analysis",
                                 "all-phenotypes-two-years",
                                 paste0("Annotation_tables_", names(df1)[2],
                                        "_", PCdf1$NumPCs, "_PCs", ".rds")))

  # Save a QQ plot
  ggsave(snp_qq(gwas[which(!is.na(gwas$score)),]),
         filename = paste0("QQplot_", names(df1)[2], "_", PCdf1$NumPCs,
                           "_PCs.png"), width = 4,
         height = 4, path = file.path("/", "home", "alice", "Github",
                                      "pvdiv-phenology-gxe", "analysis",
                                      "all-phenotypes-two-years"))

  # Save a Manhattan plot with Bonferroni
  ggsave(snp_manhattan(gwas, infos.chr = CHRN$CHRN,
                       infos.pos = snp$map$physical.pos, npoints = 1E06) +
           geom_hline(yintercept = bonferroni, lty = 2),
         filename = paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs,
                           "_PCs_bonferroni.png"), width = 12,
         height = 4, path = file.path("/", "home", "alice", "Github",
                                      "pvdiv-phenology-gxe", "analysis",
                                      "all-phenotypes-two-years"))
}
