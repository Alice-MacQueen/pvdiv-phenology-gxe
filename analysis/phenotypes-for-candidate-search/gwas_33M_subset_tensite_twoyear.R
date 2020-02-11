library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)

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
load("~/Github/pvdiv-phenology-gxe/data/phenotypes_2019.rda")

setwd("~/Github/pvdiv-phenology-gxe/analysis/all-phenotypes-two-years/")

# Prepare phenotype data frame.
df <- plants %>%
  enframe(value = "PLANT_ID") %>%
  left_join(phenotypes_2018, by = "PLANT_ID") %>%
  left_join(phenotypes_wide_2019, by = "PLANT_ID") %>%
  dplyr::select(-name)



# Look for the number of PCs to correct for population structure
# that makes lambda_GC closest to 1.05.
lambdagc <- pvdiv_lambda_GC(df = df, type = "linear", snp = snp,
                            covar = svd, ncores = NCORES,
                            npcs = c(0:8), saveoutput = TRUE)
## for CHLR, k = 3; for DAM, k = 2 and k = 5
# ------------------------------------

txdb <- loadDb(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

anno_info <-
  read_delim(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                              "Pvirgatum_516_v5.1.annotation_info.txt"),
             col_names = TRUE, delim = "\t")
gene_anno_info <- tbl_df(anno_info) %>%
  distinct(locusName, .keep_all = TRUE)

PCdf <- read_csv("~/Github/pvdiv-phenology-gxe/analysis/phenotypes-for-candidate-search/Best_Lambda_GC_KING_DAM_to_KING_CHLR_JASON_138.y_Phenotypes_0_to_8_PCs.csv")

for(i in 2:nrow(PCdf)){
  PCdf1 <- PCdf[i,]
  df1 <- df %>%
    dplyr::select(PLANT_ID, PCdf1$trait)

  gwas <- pvdiv_gwas(df = df1, type = "linear", snp = snp, covar = svd,
                     ncores = NCORES, npcs = PCdf1$NumPCs)
  # Save a Manhattan plot
  ggsave(snp_manhattan(gwas, infos.chr = CHRN$CHRN,
                       infos.pos = snp$map$physical.pos, npoints = 2E06) +
           geom_hline(yintercept = bonferroni, lty = 2),
         filename = paste0("Manhattan_", names(df1)[2], "_", PCdf1$NumPCs,
                           "_PCs.png"), width = 12,
         height = 4, path = file.path("/", "home", "alice", "Github",
                                      "pvdiv-phenology-gxe", "analysis",
                                      "phenotypes-for-candidate-search"))
}
  # Save a QQ plot
  ggsave(snp_qq(gwas[which(!is.na(gwas$score)),]),
         filename = paste0("QQplot_", names(df)[i+1], "_", PCdf1$NumPCs,
                           "_PCs.png"), width = 4,
         height = 4, path = file.path("/", "home", "alice", "Github",
                                      "pvdiv-phenology-gxe", "analysis",
                                      "phenotypes-for-candidate-search"))

  # Save annotation tables for the top associations
  anno_tables <- pvdiv_table_topsnps(df = gwas, type = "bigsnp", n = c(10,50),
                                     FDRalpha = 0.1,
                                     rangevector = c(0, 20000, 100000),
                                     markers = markers,
                                     anno_info = gene_anno_info,
                                     txdb = txdb)
  saveRDS(anno_tables, file.path("/", "home", "pubjuenger", "GWAS",
                                 paste0("Annotation_tables_", names(df)[i+1],
                                        "_", subset_names[j], "_subset_",
                                        PCdf1$NumPCs, "_PCs", ".rds")))
}
