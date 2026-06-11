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

  gen <- raw |> 
    select(
      STUDYID,
      USUBJID
    )




}




setdiff(col_by_mrinal, colnames(raw))

setdiff(col_name, colnames(raw))

sum(colnames(raw) == "ADSELVARS")

colnames(raw)[
  duplicated(col_by_mrinal, colnames(raw))
]

View(raw)


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

col_by_mrinal <- c(
'IMFL',
'PARQUAL',
'PARAMCD',
'AVAL',
'AVALC',
'AVALCAT1',
'IMEVFL',
'ANL01FL',
'ANL02FL',
'ANL03FL',
'ANL04FL',
'ANL05FL',
'ANL06FL',
'ANL07FL',
'ANL08FL',
'ANL09FL',
'ANL10FL',
'AVISIT',
'AVISITN'
) |> 
  trimws()
