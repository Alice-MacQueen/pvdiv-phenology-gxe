library(tidyverse)
library(bigsnpr)
library(switchgrassGWAS)

NCORES <- nb_cores()

snp <- snp_attach(file.path("/", "home", "alice", "Github", "pvdiv-genome",
                            "fourway_seven",
                            "Pvirgatum_V5_GWAS_630g_33M_fourway_seven.rds"))

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
                         "fourway_seven", "SVD_20PCs_630g_33M_fourway_seven.rds"))
# load("~/Github/pvdiv-phenology-gxe/data/phenotypes_2018.rda")
load("~/Github/pvdiv-phenology-gxe/data/phenotypes_2019.rda")

setwd("~/Github/pvdiv-phenology-gxe/FL50/gwas_33M/")

# Prepare phenotype data frame.
df <- plants %>%
  enframe(value = "PLANT_ID") %>%
  left_join(phenotypes_wide_2019, by = "PLANT_ID") %>%
  dplyr::select(PLANT_ID, ends_with("FL50")) %>%
  dplyr::select(-FRMI_FL50, -KING_FL50, -OVTN_FL50)

# Look for the number of PCs to correct for population structure
# that makes lambda_GC closest to 1.05.
lambdagc <- pvdiv_lambda_GC(df = df, type = "linear", snp = snp,
                            covar = svd, ncores = NCORES,
                            npcs = c(0:7), saveoutput = TRUE)
## for FL50, k = 4 and/or k = 6
# ------------------------------------

txdb <- loadDb(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

anno_info <-
  read_delim(file = file.path("/", "home", "alice", "Github", "pvdiv-genome",
                              "Pvirgatum_516_v5.1.annotation_info.txt"),
             col_names = TRUE, delim = "\t")
gene_anno_info <- tbl_df(anno_info) %>%
  distinct(locusName, .keep_all = TRUE)

PCdf <- read_csv("~/Github/pvdiv-phenology-gxe/analysis/FL50/gwas_33M/Best_Lambda_GC_BRKG_FL50_to_TMPL_FL50_Phenotypes_4_PCs.csv")

for(i in 1:nrow(PCdf)){
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
                           "_PCs.png"), width = 10,
         height = 4, path = file.path("/", "home", "alice", "Github",
                                      "pvdiv-phenology-gxe", "analysis",
                                      "FL50", "gwas_33M"))
}

# Run bigsnp2mashr

markers20per <- snp_clumping(G, infos.chr = CHRN$CHRN, ncores = NCORES)
markers2 <- markers[markers20per,] # Get unlinked markers for mash
gwas_rds <- pvdiv_results_in_folder(path = ".", pattern = "GWAS_object_") # 8 phenotypes from 4 sites - 32 phenotypes total with > 200 individuals phenotyped and sequenced.

phenotype_names <- str_sub(gwas_rds, start = 13, end = -5)

numSNPs <- ceiling(1000000/length(phenotype_names)^2)
numSNPs <- 20000 # the above gave 20409, so just rounded this to 20000.

mash_input <- pvdiv_bigsnp2mashr(path = ".", gwas_rds = gwas_rds,
                                 phenotypes = phenotype_names, numSNPs = numSNPs,
                                 markers = markers, markers2 = markers2, G = G,
                                 model = "linear", saveoutput = TRUE)

mash_output <- mash_standard_run(path = ".", list_input = mash_input,
                                 numSNPs = numSNPs, saveoutput = TRUE)

anno_tables <- pvdiv_table_topsnps(df = mash_output, type = "mash",
                                   n = c(10,500), FDRalpha = NA,
                                   rangevector = c(0, 20000, 100000),
                                   markers = markers,
                                   anno_info = gene_anno_info,
                                   txdb = txdb)
saveRDS(anno_tables, file.path("~", "Github", "pvdiv-phenology-gxe",
                               "large-files", "FL50",
                               paste0("Annotation_tables_",
                                      "Flowering_date_50per_mash_",
                                      "7_sites_4_PCs.rds")))

