
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvdiv-phenology-gxe

<!-- badges: start -->
<!-- badges: end -->

pvdiv-phenology-gxe contains the non-genomic data and code used to
analyze the genetics of GxE in the phenology data collected for the
switchgrass diversity panel in 2019, including green-up dates and
flowering dates in 2019. This repository will be private until the
associated manuscript is published.

In general, if you have any questions about how to locate an analysis
folder or code used for an analysis, look at this README first. Then, if
you still have questions, feel free to open a pull request with your
issue and I will be happy to help.

## Repository Files:

### In data-raw/:

    *metadata.R creates the metadata rda file used for subsequent analyses.
    *phenotypes_2018.R creates the phenotypic data rda file used for analyses of 2018 data
    *phenotypes_2019.Rmd creates the phenotypic data rda file used for analyses of 2019 data. These plants differ from those in phenotypes_2018 - many plants that died were replaced by Blackwell cultivar clones, and any phenotypes measured on these clones were not included in downstream analyses.
    *Taxa.R creates the Taxa rda file used to tie sequence data to the 785 sequenced PLANT_IDs.

### In data/:

    * Contains datasets as ".rds" files that are used in analyses (to produce results in the analysis/ folder). These files can be read into R using the `readRDS()` function.

### In analysis/:

    * Data objects and draft figures created as part of the analysis of heritability (including variance components analysis, done in ASReml), genome-wide association (on phenotypic BLUPs, done using the switchgrassGWAS R package), and mash (on results from the univariate GWAS, done using the switchgrassGWAS R package). Generally, the code to replicate these analyses can be found in the R/ directory.

### In R/:

    * Contains the majority of the code used in the manuscript (some additional code used to prepare data is in data-raw/.) 
    * Files that start with "Plots" were used to make manuscript figures
    * Files that start with "nohup" were used to run mash
    * Files that start with "Analysis" were used to generate files in the analysis/ folder. I tried to name these descriptively but there are a lot of them; please open a pull request if you want help locating specific code.

### In phenology-fourway/:

    * Code and data objects needed to replicate the genetic mapping of flowering-related QTL in the pseudo-F2 cross. 

### In manuscript/:

    * Contains completed drafts of a manuscript on the genetics of GxE in phenology in the switchgrass diversity panel. 
