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
gen_adishum <- function(seed = 123) {
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


# run function ------------------------------
# generate dataset
adishum <- gen_adishum()
