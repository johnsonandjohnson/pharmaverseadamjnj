# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# For `HYPSCAT` derivation. Downloaded from:
# https://www.fda.gov/drugs/development-resources/
# office-new-drugs-custom-medical-queries-ocmqs
terms <-
  file.path("data-raw", "OCMQs_v3.0.xlsm") |>
  readxl::read_excel(sheet = "Hypersensitivity") |>
  
  # For ease of operation
  dplyr::mutate(
    AEDECOD = toupper(Term),
    HYPSCAT = substr(`Algorithmic Category`, 1, 1)
  ) |>
  dplyr::select(AEDECOD, HYPSCAT)

# Generate ADAGOCMQ dataset
gen_adagocmq <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get source data
  raw <- adaeocmq
  
  # Main logic
  gen <- raw |> 
    # Derive `HYPSCAT`
    dplyr::left_join(terms)
  
  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    HYPSCAT = "Hypersensitivity Category"
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
adagocmq <- gen_adagocmq()
