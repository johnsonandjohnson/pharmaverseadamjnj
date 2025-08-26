# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADAGOCMQ dataset
gen_adagocmq <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get source data
  raw <- adaeocmq
  
  # Main logic
  
  # Define additional labels for new variables not in source dataset
  additional_labels <- list()
  
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
