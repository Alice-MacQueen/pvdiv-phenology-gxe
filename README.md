
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvdiv-phenology-gxe

<!-- badges: start -->

<!-- badges: end -->

The goal of pvdiv-phenology-gxe is to analyze the genetics of GxE in the
phenology data collected for the switchgrass diversity panel in 2018 and
2019, including tiller count in 2018, survival to 2019, and greenup to
tiller count & biomass in 2019. This repository will be private until
the associated dataset is published.

As of December 2019, this analysis is on data Jason Bonette compiled
2019-12-25. I create downstream data objects from data collected on this
date. At the beginning of every month, and whenever I hear from him, I
will check his dataset for any changes.

R Files:

In data-raw/:

    *metadata.R creates the metadata rda file used for subsequent analyses.
    *phenotypes_2018.R creates the phenotypic data rda file used for analyses of 2018 data
    *phenotypes_2019.Rmd creates the phenotypic data rda file used for analyses of 2019 data. These plants differ from those in phenotypes_2018 - many plants that died were replaced by Blackwell cultivar clones, and any phenotypes measured on these clones were not included in downstream analyses.
    *Taxa.R creates the Taxa rda file used to tie sequence data to the 785 sequenced PLANT_IDs.
