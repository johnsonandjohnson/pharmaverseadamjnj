# Generate ADPC dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADPC dataset
gen_adpc <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get source data
  raw <- pharmaverseadam::adpc

  gen <- dplyr::mutate(
    raw,
    CRIT1 = dplyr::if_else(grepl("<", PCSTRESC), "<LLOQ", NA),
    CRIT1FL = dplyr::if_else(grepl("<", PCSTRESC), "Y", "N")
  )
  
  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    CRIT1 = "Concentrations Below the Lowest Quantification Level Criterion",
    CRIT1FL = "Concentrations Below the Lowest Quantification Level Flag"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  return(gen)
}

# Generate the dataset
adpc <- gen_adpc()
