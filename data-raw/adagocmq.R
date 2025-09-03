# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))
source(file.path("data-raw", "adagocmq-helpers.R"))


# Generate ADAGOCMQ dataset
gen_adagocmq <- function() {
  # Get source data
  raw_adaeocmq <- adaeocmq
  raw_adlb <- adlb
  raw_adcm <- adcm
  
  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    HYPSCAT = "Hypersensitivity Category",
    ATERMN = "Analysis Term",
    ATERM = "Analysis Term (N)"
  )
  
  
  # Create records related to ADAEOCMQ ----------------------------------------

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
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypersensitivity" & HYPSCAT == "A") |>
      dplyr::mutate(ATERMN = 11),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypersensitivity" & HYPSCAT == "B") |>
      dplyr::mutate(ATERMN = 121),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypersensitivity" & HYPSCAT == "C") |>
      dplyr::mutate(ATERMN = 122),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypersensitivity" & HYPSCAT == "D") |>
      dplyr::mutate(ATERMN = 131),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hyperglycemia" & OCMQCLSS == "Narrow") |>
      dplyr::mutate(ATERMN = 21),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypoglycemia" & OCMQCLSS == "Narrow") |>
      dplyr::mutate(ATERMN = 31),
    
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Hypoglycemia" & OCMQCLSS == "Broad") |>
      dplyr::filter(AEDECOD %in% toupper(c(
        "Accident", "Anxiety", "Asthenia", "Cold sweat", "Coma",
        "Confusional state", "Fall", "Fatigue", "Hunger", "Hyperhidrosis",
        "Irritability", "Loss of consciousness", "Palpitations",
        "Road traffic accident", "Seizure", "Tremor", "Dysarthria",
        "Balance disorder", "Coordination abnormal", "Headache",
        "Vision blurred", "Visual impairment"
      ))) |>
      dplyr::mutate(ATERMN = 331)
  )
  
  
  # Create records related to ADLB --------------------------------------------

  # Derive `ATERMN`
  records_adlb <- dplyr::bind_rows(
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y") |>
      dplyr::filter(grepl("mmol/L", PARAM) & AVAL >= 6.993) |>
      dplyr::mutate(ATERMN = 22),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y") |>
      dplyr::filter(grepl("mg/dL", PARAM) & AVAL >= 126) |>
      dplyr::mutate(ATERMN = 22),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & grepl("mmol/L", PARAM) & AVAL > 9.99) |>
      dplyr::mutate(ATERMN = 231),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & grepl("mg/dL", PARAM) & AVAL > 180) |>
      dplyr::mutate(ATERMN = 231),

    raw_adlb |>
      dplyr::filter(PARAMCD == "HBA1C" & AVAL >= 6.5 & ADT >= TRTSDT) |>
      dplyr::mutate(ATERMN = 25),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "HBA1C" & AVAL >= 5.7 & CHG >= 0.3) |>
      dplyr::mutate(ATERMN = 26),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y") |>
      dplyr::filter(grepl("mmol/L", PARAM) & CHG >= 1.11 & AVAL > 5.55) |>
      dplyr::mutate(ATERMN = 27),
    
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA" & LBFAST == "Y") |>
      dplyr::filter(grepl("mg/dL", PARAM) & CHG >= 20 & AVAL > 100) |>
      dplyr::mutate(ATERMN = 27),

    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA") |>
      dplyr::filter(grepl("mmol/L", PARAM) & AVAL < 3.0) |>
      dplyr::mutate(ATERMN = 32),

    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA") |>
      dplyr::filter(grepl("mg/dL", PARAM) & AVAL < 54) |>
      dplyr::mutate(ATERMN = 32)
  )
  
  
  # Create records related to ADCM --------------------------------------------
  
  # Derive `ATERMN = 24`
  records_adcm <- raw_adcm |>
    dplyr::filter(ASTDT >= TRTSDT) |>
    dplyr::filter(grepl(
      paste(
        "diab", "mellitus", "hyperglyc", "glucose", "dibet", "dieb",
        sep = "|"
      ),
      CMINDC,
      ignore.case = TRUE
    )) |>
    dplyr::filter(!grepl(
      paste(
        "prophyla", "prevent", "insipidus", "hyperglycerid",
        "low blood glucose", "low glucose", "low blood sugar", "low sugar",
        "low afternoon blood glucose", "low morning blood glucose",
        sep = "|"
      ),
      CMINDC,
      ignore.case = TRUE
    )) |>
    dplyr::filter(grepl(
      paste(
        "gliptin", "glutide", "diabet", "glitaz", "glucose lowering",
        "glucosidas", "dipeptidyl", "sulfonyl", "DPP", "guanide", "GLP",
        "glucagon-like", "metform", "gliflozin", "insulin", "sodium-glucose",
        "SGLT", "thiazolid",
        sep = "|"
      ),
      CMCLAS,
      ignore.case = TRUE
    )) |>
    dplyr::filter(!grepl("sex hormone", CMCLAS, ignore.case = TRUE)) |>
    dplyr::mutate(ATERMN = 24)
  
  
  # Combine records -----------------------------------------------------------
  
  gen <- dplyr::bind_rows(records_adaeocmq, records_adlb, records_adcm) |>
    # Derive new records based on `ATERMN`
    derive_atermn_1x(c(121, 122), 12) |>
    derive_atermn_1x(c(121, 131), 13) |>
    derive_atermn_1x(c(122, 131), 14) |>
    derive_atermn_23() |>
    
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
    ))
  
  return(gen)
}


# Generate the dataset
adagocmq <- gen_adagocmq()
