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
    CRIT1 = "<LLOQ",
    CRIT1FL = dplyr::if_else(grepl("<", PCSTRESC), "Y", "N"),
    PKRFL = dplyr::case_when(
      # This condition checks if a record should be flagged as "Y" (Yes) for PK evaluable
      # It requires ALL of the following to be true:
      # 1. None of the CRIT flags (CRIT1FL, CRIT2FL, etc.) are "Y" - meaning no criterion violations
      # 2. PCDTC (date/time) is not missing (NA)
      # 3. PCDTC matches the ISO 8601 datetime format (YYYY-MM-DDTHH:MM) - ensuring valid timestamp
      !dplyr::if_any(dplyr::matches("^CRIT\\d+FL$"), ~ . == "Y") &
        !is.na(PCDTC) & grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}", PCDTC) ~ "Y",
      # This condition checks if a record should be flagged as "N" (No) for PK evaluable
      # It requires ANY of the following to be true:
      # 1. At least one of the CRIT flags (CRIT1FL, CRIT2FL, etc.) is "Y" - meaning there is a criterion violation
      # 2. PCDTC (date/time) is missing (NA)
      # 3. PCDTC does NOT match the ISO 8601 datetime format (YYYY-MM-DDTHH:MM) - meaning invalid timestamp
      dplyr::if_any(dplyr::matches("^CRIT\\d+FL$"), ~ . == "Y") |
        is.na(PCDTC) | !grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}", PCDTC) ~ "N",
      .default = NA_character_
    )
  )

  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    CRIT1 = "Lowest Quantification Level Criterion",
    CRIT1FL = "Lowest Quantification Level Flag",
    PKRFL = "PK Analysis Record-Level Flag"
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
