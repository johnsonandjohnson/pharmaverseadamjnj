# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))
source(file.path("data-raw", "adagocmq-helpers.R"))


# Generate ADAGOCMQ dataset
gen_adagocmq <- function() {
  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    HYPSCAT = "Hypersensitivity Category",
    ATERMN = "Analysis Term",
    ATERM = "Analysis Term (N)"
  )
  
  # Create records related to ADAEOCMQ ----------------------------------------

  # Get source data
  raw_adaeocmq <- adaeocmq
  
  # For `HYPSCAT` derivation. Downloaded from:
  # https://www.fda.gov/drugs/development-resources/
  # office-new-drugs-custom-medical-queries-ocmqs
  terms <-
    file.path("data-raw", "OCMQs_v3.0.xlsm") |>
    readxl::read_excel(sheet = "Hypersensitivity") |>
    dplyr::mutate(
      AEDECOD = toupper(Term),
      HYPSCAT = substr(`Algorithmic Category`, 1, 1)
    ) |>
    dplyr::select(AEDECOD, HYPSCAT)
  
  # Derive `HYPSCAT`
  raw_adaeocmq <- raw_adaeocmq |>
    dplyr::left_join(terms)
  
  # Derive `ATERMN`
  records_adaeocmq <- dplyr::bind_rows(
    derive_atermn(
      raw_adaeocmq,
      rule = 'OCMQNAM == "Hypersensitivity" & HYPSCAT == "A"',
      level = 11
    ),
    derive_atermn(
      raw_adaeocmq,
      rule = 'OCMQNAM == "Hypersensitivity" & HYPSCAT == "B"',
      level = 121
    ),
    derive_atermn(
      raw_adaeocmq,
      rule = 'OCMQNAM == "Hypersensitivity" & HYPSCAT == "C"',
      level = 122
    ),
    derive_atermn(
      raw_adaeocmq,
      rule = 'OCMQNAM == "Hypersensitivity" & HYPSCAT == "D"',
      level = 131
    ),
    derive_atermn(
      raw_adaeocmq,
      rule = 'OCMQNAM == "Hyperglycemia" & OCMQCLSS == "Narrow"',
      level = 21
    )
  )
  
  records_adaeocmq <- records_adaeocmq |>
    # Derive new records based on `ATERMN`
    derive_atermn_1x(c(121, 122), 12) |>
    derive_atermn_1x(c(121, 131), 13) |>
    derive_atermn_1x(c(122, 131), 14) |>
    
    # Derive `ATERM`
    dplyr::mutate(ATERM = dplyr::case_when(
      ATERMN == 11 ~ "Any hypersensitivity OCMQ narrow term",
      ATERMN == 121 ~ "Respiratory",
      ATERMN == 122 ~ "Skin Reaction",
      ATERMN == 12 ~ "Respiratory + Skin Reaction",
      ATERMN == 131 ~ "Systemic Reaction",
      ATERMN == 13 ~ "Respiratory + Systemic Reaction",
      ATERMN == 14 ~ "Skin + Systemic Reaction",
      ATERMN == 21 ~ "Any Hyperglycemia OCMQ Narrow Term"
    )) |>
    
    # Handle NA values and convert characters to factors
    df_na(char_as_factor = TRUE) |>
    # Restore labels
    restore_labels(raw_adaeocmq, additional_labels)
  
  
  # Create records related to ADLB --------------------------------------------
  
  # Get source data
  raw_adlb <- adlb
  
  # Derive `ATERMN`
  records_adlb <- dplyr::bind_rows(
    derive_atermn(
      raw_adlb,
      rule = paste(
        'PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y" &',
        'grepl("mmol/L", PARAM) & AVAL >= 6.993'
      ),
      level = 22
    ),
    derive_atermn(
      raw_adlb,
      rule = paste(
        'PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y" &',
        'grepl("mg/dL", PARAM) & AVAL >= 126'
      ),
      level = 22
    ),
    derive_atermn(
      raw_adlb,
      rule = 'PARAMCD == "GLUC" & grepl("mmol/L", PARAM) & AVAL > 9.99',
      level = 231
    ),
    derive_atermn(
      raw_adlb,
      rule = 'PARAMCD == "GLUC" & grepl("mg/dL", PARAM) & AVAL > 180',
      level = 231
    ),
    derive_atermn(
      raw_adlb,
      rule = 'PARAMCD == "HBA1C" & AVAL >= 6.5 & ADT >= TRTSDT',
      level = 25
    ),
    derive_atermn(
      raw_adlb,
      rule = 'PARAMCD == "HBA1C" & AVAL >= 5.7 & CHG >= 0.3',
      level = 26
    ),
    derive_atermn(
      raw_adlb,
      rule = paste(
        'PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y" &',
        'grepl("mmol/L", PARAM) & CHG >= 1.11 & AVAL > 5.55'
      ),
      level = 27
    ),
    derive_atermn(
      raw_adlb,
      rule = paste(
        'PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y" &',
        'grepl("mg/dL", PARAM) & CHG >= 20 & AVAL > 100'
      ),
      level = 27
    )
  )
  
  records_adlb <- records_adlb |>
    # Derive new records based on `ATERMN`
    derive_atermn_23()
  
  
  return(gen)
}


# Generate the dataset
adagocmq <- gen_adagocmq()
