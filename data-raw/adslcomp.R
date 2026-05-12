# Generate ADSLCOMP dataset

# Load necessary libraries
library(dplyr)
library(junco)
library(tibble)


# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADSLCOMP dataset
gen_adslcomp <- function(seed = 123) {
  # Source related datasets
  source(file.path("data-raw", "adsl.R"))

  # Generate ADSL with "Std Of Care" arm
  gen <- make_fake_adsl(adsl)

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels using raw and additional_labels
  gen <- restore_labels(
    df = gen,
    orig_df = adsl, # Use original adsl
  )

  return(gen)
}

adslcomp <- gen_adslcomp()
