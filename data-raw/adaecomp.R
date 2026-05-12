# Generate ADAECOMP dataset

# Load necessary libraries
library(dplyr)
library(junco)
library(tibble)


# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADAECOMP dataset
gen_adaecomp <- function(seed = 123) {
  # Source related datasets
  source(file.path("data-raw", "adsl.R"))
  source(file.path("data-raw", "adae.R"))

  # Generate ADSL with "Std Of Care" arm
  adsl <- make_fake_adsl(adsl)

  # Generate ADAE with borrowed events
  gen <- adae |>
    filter(TRTEMFL == "Y") |>
    borrow_aes(adsl)

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels using raw and additional_labels
  gen <- restore_labels(
    df = gen,
    orig_df = adae, # Use original adae
  )

  return(gen)
}

adaecomp <- gen_adaecomp()
