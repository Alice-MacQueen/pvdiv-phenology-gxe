
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvdiv-phenology-gxe

<!-- badges: start -->
<!-- badges: end -->

pvdiv-phenology-gxe contains the non-genomic data and code used to
analyze the genetics of GxE in the phenology data collected for the
switchgrass diversity panel in 2019, including green-up dates (the
timing of vegetative growth) and flowering dates (the timing of
reproductive growth) in 2019.

The **most recent submission** is on bioRxiv:
<https://doi.org/10.1101/2021.08.19.456975>

**Code** used in the 2024 submission is located here:
<https://github.com/Alice-MacQueen/pvdiv-phenology-gxe/blob/main/R/Analysis_v1.2_mash_using_greedy_covar.qmd>.

You can find the **genomic data** at:
<https://doi.org/10.18738/T8/A604BU>.

In general, if you have any questions about how to locate an analysis
folder or code used for an analysis, look at this README first. Then, if
you still have questions, feel free to open a pull request with your
issue or contact me directly at alice.macqueen\[at\]gmail.com and I will
be happy to help.

## Repository Files:

### In data-raw/:

    - metadata.R creates the metadata rda file used for subsequent analyses.
    - phenotypes_2018.R creates the phenotypic data rda file used for analyses of 2018 data
    - phenotypes_2019.Rmd creates the phenotypic data rda file used for analyses of 2019 data. These plants differ from those in phenotypes_2018 - many plants that died were replaced by Blackwell cultivar clones, and any phenotypes measured on these clones were not included in downstream analyses.
    - Taxa.R creates the Taxa rda file used to tie sequence data to the 785 sequenced PLANT_IDs.

### In data/:

    - Contains datasets as ".rds" files that are used in analyses (to produce results in the analysis/ folder). These files can be read into R using the `readRDS()` function.

### In analysis/:

Data objects and draft figures created as part of different discrete
analyses, including:

    - heritability/ contains the analysis of heritability (including variance components analysis, done in ASReml);
    - gwas/BLUPs/  contains genome-wide association results (on phenotypic BLUPs, done using the switchgrassGWAS R package);
    - hypothesis_matrices/ contains unpublished analyses of correlation data between GxWeather matrices, used as exploratory plots for a revision;
    - mash/ contains mash models (on results from the univariate GWAS, done using the switchgrassGWAS R package). 
    - mash_greedy_algorithm/mash_h2_diag_covar/ contains runs of the greedy mash algorithm used to select covariance matrices for inclusion in the mash model for the current revision. mash_greedy_algorithm/ also contains other runs of the greedy mash algorithm used as exploratory analyses for suitability for the revision.

Generally, the code to replicate these analyses can be found in the R/
directory.

### In R/:

R/ contains the majority of the code used in the manuscript (some
additional code used to prepare data is in data-raw/.)

    - The file titled "Analysis_v1.2_mash_using_greedy_covar.qmd" contains all of the code used to run the analyses and generate the figures for the 2024 revision (in manuscript/revision2/). 
    - THe file titled "SI_Appendix_Datasets.Rmd" contains the code used to generate the SI Appendix datasets and the SI Appendix figures for the 2024 revision.
    - The file titled "UT_Dataverse_files.Rmd" was used to generate the SNP datasets uploaded to the UT Dataverse for replication of this analysis. 
    - Files that start with "Plots" were used to make manuscript figures for a 2021 submission
    - Files that start with "nohup" were used to run mash for a 2021 submission
    - Files that start with "Analysis" were used to generate files in the analysis/ folder. I tried to name these descriptively but there are a lot of them; please open a pull request if you want help locating specific code from earlier versions of this analysis. The most recent version is in the file titled "Analysis_v1.2_mash_using_greedy_covar.qmd".

### In phenology-fourway/:

    - Code and data objects needed to replicate the genetic mapping of flowering-related QTL in the pseudo-F2 cross. 

### In manuscript/:

This folder contains completed drafts of a manuscript on the genetics of
GxE in phenology in the switchgrass diversity panel.

    - revision2/ contains the most recent, 2024 revision. The main text is called revision2.pdf. There are also folders containing pdfs of the response to reviews (in response_to_reviews/), the SI Appendix data (in SI_Appendix) and images used in the main manuscript (in images/).
    - Coauthors/ contains the manuscript from a submission in 2021.
