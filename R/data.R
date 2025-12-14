#' Small OTU abundance dataset
#'
#' A dataset containing relatives abundances of 60 OTUs for 60 human stool samples.
#' This is a subset of the data provided in `extdata/otu_large.csv`, which was
#' used in [Topçuoğlu _et al._ 2020](https://journals.asm.org/doi/10.1128/mbio.00434-20).
#'
#' @format A data frame with 60 rows and 61 variables.
#' The `dx` column is the diagnosis: healthy or cancerous (colorectal).
#' All other columns are OTU relative abundances.
"pop_data_us_23"

#' Mini OTU abundance dataset
#'
#' A dataset containing relatives abundances of OTUs for human stool samples
#' with a binary outcome, `dx`.
#' This is a subset of `otu_small`.
#'
#' @format A data frame
#' The `dx` column is the diagnosis: healthy or cancerous (colorectal).
#' All other columns are OTU relative abundances.
"sim_sp_ef"

#' Mini OTU abundance dataset - preprocessed
#'
#' This is the result of running `preprocess_data("otu_mini_bin")`
"temp_data_us_23"

#' Mini OTU abundance dataset with 3 categorical variables
#'
#' A dataset containing relatives abundances of OTUs for human stool samples
#'
#' @format A data frame
#' The `dx` column is the colorectal cancer diagnosis: adenoma, carcinoma, normal.
#' All other columns are OTU relative abundances.
"us_county_23"

