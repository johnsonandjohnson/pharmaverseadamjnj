# generate adishum dataset
# setup section ------------------------------

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# core function ------------------------------

# Generate adishum dataset
get_adishum <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # get source data
  raw <- pharmaverseadam::adis_vaccine

  gen <- raw() |> 
    select(
      STUDYID,
      USUBJID
    )




}


names(raw_df)
col_name

col_name <- c(
  'STUDYID',
  'USUBJID',
  'ASEQ',
  'ADT',
  'ATM',
  'ADTM',
  'ADY',
  'APERDY',
  'APHDY',
  'AVISIT',
  'AVISITN',
  'APERIOD',
  'APERIODC',
  'APERSDT',
  'APERSDTM',
  'APEREDT',
  'APEREDTM',
  'APHASE',
  'APHASEN',
  'PHSDT'
)
