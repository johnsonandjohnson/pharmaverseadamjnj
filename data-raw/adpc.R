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
  
  gen <- dplyr::mutate(
    raw,
    
    PCTESTCD = as.factor(dplyr::case_match(
      PCTESTCD,
      "XAN" ~ "APA",
      .default = PCTESTCD
    )),
    
    PCTEST = as.factor(dplyr::case_match(
      PCTEST,
      "XANOMELINE" ~ "APALUTAMIDE",
      .default = PCTEST
    )),
    
    TRT01P = forcats::fct_recode(
      TRT01P,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    
    TRT01A = forcats::fct_recode(
      TRT01A,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    
    PARAMCD = as.factor(dplyr::case_match(
      PARAMCD,
      "XAN" ~ "APA",
      .default = PARAMCD
    )),
    
    PARAM = as.factor(gsub("Xanomeline", "Apalutamide", PARAM)),

    ARMCD = as.factor(dplyr::case_when(
      ARMCD == "Xan_Hi" ~ "Apa",
      ARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ARMCD
    )),
    
    ARM = as.factor(dplyr::case_when(
      ARM == "Xanomeline High Dose" ~ "Apalutamide",
      ARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ARM
    )),
    
    ACTARMCD = as.factor(dplyr::case_when(
      ACTARMCD == "Xan_Hi" ~ "Apa",
      ACTARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ACTARMCD
    )),
    
    ACTARM = as.factor(dplyr::case_when(
      ACTARM == "Xanomeline High Dose" ~ "Apalutamide",
      ACTARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ACTARM
    ))
  )
  
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
