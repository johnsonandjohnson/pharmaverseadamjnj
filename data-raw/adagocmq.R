# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))
# Source utility functions for creating new records based on complex rules
source(file.path("data-raw", "adagocmq-helpers.R"))

# Source related datasets
source(file.path("data-raw", "adaeocmq.R"))
source(file.path("data-raw", "adlb.R"))
source(file.path("data-raw", "adcm.R"))
source(file.path("data-raw", "adsl.R"))


#' Generate ADAGOCMQ Dataset
#'
#' The code below has five parts:
#'
#' 1. Create new records related to ADAEOCMQ.
#' 2. Create new records related to ADLB.
#' 3. Create new records related to ADCM.
#' 4. Copy variables from ADSL.
#' 5. Combine data and derive variables.
#'
#' The creation of new records is driven by the derivation of variable `ATERMN`.
#' The file adagocmq-helpers.R contains utility functions for deriving `ATERMN`
#' based on complex rules.
#'
#' @noRd
gen_adagocmq <- function() {
  # Variables to keep for created records
  records_variables <- c(
    "USUBJID",
    "ASTDT",
    "ASTDY",
    "ATERMN",
    "SRCVALUE",
    "SRCVAR",
    "SRCDOM",
    "SRCSEQ"
  )

  # Create records related to ADAEOCMQ ----------------------------------------

  # Get source data
  raw_adaeocmq <- adaeocmq

  # For `HYPSCAT` derivation. Downloaded from:
  # https://www.fda.gov/drugs/development-resources/
  # office-new-drugs-custom-medical-queries-ocmqs
  terms <-
    file.path("source_data", "OCMQs_v3.0.xlsm") |>
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
      dplyr::mutate(ATERMN = 331),
    raw_adaeocmq |>
      dplyr::filter(OCMQNAM == "Muscle Injury" & OCMQCLSS == "Narrow") |>
      dplyr::mutate(ATERMN = 41),
    raw_adaeocmq |>
      dplyr::filter(AEDECOD == toupper("Myalgia")) |>
      dplyr::mutate(ATERMN = 441),
    raw_adaeocmq |>
      dplyr::filter(AEDECOD == toupper("Muscular Weakness")) |>
      dplyr::mutate(ATERMN = 442),
    raw_adaeocmq |>
      dplyr::filter(
        AEDECOD %in% toupper(c("Myoglobin Urine Present", "Chromaturia"))
      ) |>
      dplyr::mutate(ATERMN = 443)
  ) |>
    dplyr::mutate(
      SRCVALUE = AEDECOD,
      SRCVAR = "AEDECOD",
      SRCDOM = "ADAEOCMQ",
      SRCSEQ = AESEQ
    ) |>
    dplyr::select(dplyr::all_of(c(records_variables, "HYPSCAT", "AEDECOD")))


  # Create records related to ADLB --------------------------------------------

  # Get source data
  raw_adlb <- adlb

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
      dplyr::mutate(ATERMN = 32),
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA") |>
      dplyr::filter(grepl("mmol/L", PARAM) & AVAL < 3.885) |>
      dplyr::mutate(ATERMN = 332),
    raw_adlb |>
      dplyr::filter(PARAMCD == "GLUC" & LBSPEC == "PLASMA") |>
      dplyr::filter(grepl("mg/dL", PARAM) & AVAL < 70) |>
      dplyr::mutate(ATERMN = 332),
    raw_adlb |>
      dplyr::filter(PARAMCD == "MGB" & PARCAT1 == "URINALYSIS") |>
      dplyr::filter(AVAL > ANRHI) |>
      dplyr::mutate(ATERMN = 42),
    create_atermn_43(raw_adlb)
  ) |>
    dplyr::mutate(
      ASTDT = ADT,
      ASTDY = ADY,
      SRCVALUE = as.character(AVAL),
      SRCVAR = "AVAL",
      SRCDOM = "ADLB",
      SRCSEQ = ASEQ
    ) |>
    dplyr::select(dplyr::all_of(c(records_variables, "PARAM")))


  # Create records related to ADCM --------------------------------------------

  # Get source data
  raw_adcm <- adcm

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
    dplyr::mutate(ATERMN = 24) |>
    dplyr::mutate(
      SRCVALUE = CMDECOD,
      SRCVAR = "CMDECOD",
      SRCDOM = "ADCM",
      SRCSEQ = CMSEQ
    ) |>
    dplyr::select(dplyr::all_of(records_variables))


  # Copy variables from ADSL --------------------------------------------------

  # Get source data
  raw_adsl <- adsl

  variables_adsl <- raw_adsl |>
    dplyr::select(
      STUDYID,
      USUBJID,
      TRT01P,
      TRT01PN,
      TRT01A,
      TRT01AN,
      AGE,
      AGEU,
      SEX,
      RACE,
      COUNTRY,
      RANDFL,
      SAFFL,
      SITEID,
      SUBJID
    )


  # Combine records -----------------------------------------------------------

  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    ASTDT = "Analysis Start Date",
    ASTDY = "Analysis Start Relative Day",
    ACAT1 = "Analysis Category 1",
    ACAT1N = "Analysis Category 1 (N)",
    ATERM = "Analysis Term (N)",
    ATERMN = "Analysis Term",
    HYPSCAT = "Hypersensitivity Category",
    SRCVALUE = "Source Value",
    SRCVAR = "Source Variable",
    SRCDOM = "Source Data",
    SRCSEQ = "Source Sequence Number",
    ANL01FL = "Analysis Flag 01"
  )

  gen <- dplyr::bind_rows(records_adaeocmq, records_adlb, records_adcm) |>
    # Derive new records based on `ATERMN`
    derive_atermn_1x(c(121, 122), 12) |>
    derive_atermn_1x(c(121, 131), 13) |>
    derive_atermn_1x(c(122, 131), 14) |>
    derive_atermn_23() |>
    derive_atermn_1x(c(331, 332), 33) |>
    derive_atermn_34() |>
    derive_atermn_44() |>
    # Derive `ATERM`
    dplyr::mutate(ATERM = dplyr::case_when(
      ATERMN == 11 ~ "Any hypersensitivity OCMQ narrow term",
      ATERMN == 121 ~ "Respiratory",
      ATERMN == 122 ~ "Skin Reaction",
      ATERMN == 12 ~ "Respiratory + Skin Reaction",
      ATERMN == 131 ~ "Systemic Reaction",
      ATERMN == 13 ~ "Respiratory + Systemic Reaction",
      ATERMN == 14 ~ "Skin + Systemic Reaction",
      ATERMN == 21 ~ "Any Hyperglycemia OCMQ Narrow Term",
      ATERMN == 22 & grepl("mg/dL", PARAM) ~
        "Fasting Plasma Glucose >= 126 mg/dL",
      ATERMN == 22 & grepl("mmol/L", PARAM) ~
        "Fasting Plasma Glucose >= 6.993 mmol/L",
      ATERMN == 231 & grepl("mg/dL", PARAM) ~ "Plasma Glucoses > 180 mg/dL",
      ATERMN == 231 & grepl("mmol/L", PARAM) ~ "Plasma Glucoses > 9.99 mmol/L",
      ATERMN == 23 & grepl("mg/dL", PARAM) ~ ">= 2 Plasma Glucoses > 180 mg/dL",
      ATERMN == 23 & grepl("mmol/L", PARAM) ~
        ">= 2 Plasma Glucoses > 9.99 mmol/L",
      ATERMN == 24 ~ "Any New Diabetes Concomitant Medication",
      ATERMN == 25 ~ "Post Baseline HbA1c >= 6.5%",
      ATERMN == 26 ~ "HbA1c Increase >= 0.3% with Post Baseline HbA1c >= 5.7%",
      ATERMN == 27 & grepl("mg/dL", PARAM) ~ paste(
        "Change from Baseline Fasting Plasma Glucose >= 20 mg/dL",
        "with Post Baseline Fasting Plasma Glucose > 100 mg/dL"
      ),
      ATERMN == 27 & grepl("mmol/L", PARAM) ~ paste(
        "Change from Baseline Fasting Plasma Glucose >= 1.11 mmol/L",
        "with Post Baseline Fasting Plasma Glucose > 5.55 mmol/L"
      ),
      ATERMN == 31 ~ "Any Hypoglycemia OCMQ narrow term",
      ATERMN == 32 & grepl("mg/dL", PARAM) ~ "Plasma Glucose < 54 mg/dL",
      ATERMN == 32 & grepl("mmol/L", PARAM) ~ "Plasma Glucose < 3.0 mmol/L",
      ATERMN == 331 ~ "Hypoglycemia Term",
      ATERMN == 332 & grepl("mg/dL", PARAM) ~ "Plasma Glucose < 70 mg/dL",
      ATERMN == 332 & grepl("mmol/L", PARAM) ~ "Plasma Glucose < 3.885 mmol/L",
      ATERMN == 33 & grepl("mg/dL", PARAM) ~
        "Hypoglycemia Term + Plasma Glucose < 70 mg/dL",
      ATERMN == 33 & grepl("mmol/L", PARAM) ~
        "Hypoglycemia Term + Plasma Glucose < 3.885 mmol/L",
      ATERMN == 34 & grepl("mg/dL", PARAM) ~
        ">= 2 Hypoglycemia Terms + >= 2 Episodes of Plasma Glucose < 70 mg/dL",
      ATERMN == 34 & grepl("mmol/L", PARAM) ~ paste(
        ">= 2 Hypoglycemia Terms +",
        ">= 2 Episodes of Plasma Glucose < 3.885 mmol/L"
      ),
      ATERMN == 41 ~ "Any Muscle Injury OCMQ Narrow",
      ATERMN == 42 ~ "Urine myoglobin > ULN",
      ATERMN == 43 ~ "CPK >5 x ULN",
      ATERMN == 441 ~ AEDECOD,
      ATERMN == 442 ~ AEDECOD,
      ATERMN == 443 ~ "Myoglobin Urine Present or Chromaturia",
      ATERMN == 44 ~ "Myalgia + Weakness + Chromaturia"
    )) |>
    # Derive `ACAT1N` and `ACAT1`
    dplyr::mutate(ACAT1N = substr(ATERMN, 1, 1)) |>
    dplyr::mutate(ACAT1 = dplyr::case_when(
      ACAT1N == 1 ~ "Hypersensitivity",
      ACAT1N == 2 ~ "Hyperglycemia",
      ACAT1N == 3 ~ "Hypoglycemia",
      ACAT1N == 4 ~ "Rhabdomyolysis"
    )) |>
    # Derive `ANL01FL`
    dplyr::mutate(ANL01FL = dplyr::if_else(nchar(ATERMN) == 2, "Y", NA)) |>
    # Copy and select variables
    dplyr::left_join(variables_adsl) |>
    dplyr::select(
      STUDYID,
      USUBJID,
      ASTDT,
      ASTDY,
      TRT01P,
      TRT01PN,
      TRT01A,
      TRT01AN,
      ACAT1,
      ACAT1N,
      ATERM,
      ATERMN,
      HYPSCAT,
      SRCVALUE,
      SRCVAR,
      SRCDOM,
      SRCSEQ,
      ANL01FL,
      AGE,
      AGEU,
      SEX,
      RACE,
      COUNTRY,
      RANDFL,
      SAFFL,
      SITEID,
      SUBJID
    ) |>
    # Handle NA values and convert characters to factors
    df_na(char_as_factor = TRUE) |>
    # Restore labels
    restore_labels(raw_adsl, additional_labels)

  return(gen)
}


# Generate the dataset
adagocmq <- gen_adagocmq()
