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

  # Build the same extended ADSL as adslcomp (2 non-active + 3 active arms)
  adsl_ext <- make_fake_adsl(adsl, newlev = "Std Of Care", oldlev = "Placebo")
  adsl_ext <- make_fake_adsl(adsl_ext, newlev = "Xanomeline Mid Dose", oldlev = "Xanomeline Low Dose")

  # Borrow AEs for both new arms
  gen <- adae |>
    filter(TRTEMFL == "Y") |>
    borrow_aes(adsl_ext, newlev = "Std Of Care", oldlev = "Placebo") |>
    borrow_aes(adsl_ext, newlev = "Xanomeline Mid Dose", oldlev = "Xanomeline Low Dose")

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels using raw and additional_labels
  gen <- restore_labels(
    df = gen,
    orig_df = adae,
  )

  return(gen)
}

adaecomp <- gen_adaecomp()
