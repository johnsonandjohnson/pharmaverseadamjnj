# Generate ADPC dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(forcats)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADPC dataset
gen_adpc <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get source data
  raw <- pharmaverseadam::adpc
  
  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw
  )

  return(gen)
}

# Generate the dataset
adpc <- gen_adpc()
