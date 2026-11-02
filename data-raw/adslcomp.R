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

  # Add second non-active arm ("Std Of Care") by duplicating Placebo
  gen <- make_fake_adsl(adsl, newlev = "Std Of Care", oldlev = "Placebo")

  # Add third active arm by duplicating an existing active arm
  gen <- make_fake_adsl(gen, newlev = "Xanomeline Mid Dose", oldlev = "Xanomeline Low Dose")

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels using raw and additional_labels
  gen <- restore_labels(
    df = gen,
    orig_df = adsl,
  )

  return(gen)
}

adslcomp <- gen_adslcomp()
